"""变异校验与参考左对齐归一化（纯函数，不依赖 Flask）。

支持的变异形式（均为 ACGT 大写、1-based 坐标）：
  * 单碱基替换：len(REF)==1 且 len(ALT)==1，REF != ALT
  * 带锚碱基的纯插入：len(REF)==1 且 len(ALT)>=2，且 ALT[0] == REF
  * 带锚碱基的纯缺失：len(REF)>=2 且 len(ALT)==1，且 REF[0] == ALT

归一化采用经典 left-align 流程（与 bcftools norm / vt normalize 等价）：
  1. 不断去掉两个等位基因相同的末尾碱基；
  2. 无法去尾时，把事件整体左移一位，在两个等位基因前补参考前缀碱基；
  3. 到达参考起点后，若出现空等位基因则向右补一个参考碱基重新锚定；
  4. 最后去除共同前缀（每条等位基因至少保留 1 个碱基）。
同一次插入/缺失在重复序列中无论以何种坐标书写，归一化结果都唯一。
"""

from __future__ import annotations

ACGT = frozenset("ACGT")
VALID_GTS = ("0/1", "1/1")
MAX_NAME_LEN = 128


class ValidationError(Exception):
    """输入数据不合法；message 可直接展示给用户。"""


# --------------------------------------------------------------------------- #
# 参考序列
# --------------------------------------------------------------------------- #
def parse_reference(raw: bytes | str) -> str:
    """解析参考序列上传内容，支持纯文本和单记录 FASTA。仅接受大写 ACGT。"""
    text = _decode(raw)
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValidationError("参考序列为空。")

    is_fasta = any(line.startswith(">") for line in lines)
    if is_fasta:
        headers = [line for line in lines if line.startswith(">")]
        if len(headers) != 1:
            raise ValidationError("FASTA 参考序列只能包含一条记录（仅允许一个 '>' 标题行）。")
        seq = "".join(line for line in lines if not line.startswith(">"))
    else:
        seq = "".join(lines)  # 去掉所有空白/换行

    if not seq:
        raise ValidationError("参考序列为空。")
    bad = sorted(set(seq) - ACGT)
    if bad:
        raise ValidationError(
            "参考序列只允许包含大写 A/C/G/T，发现非法字符：" + " ".join(bad[:10])
        )
    return seq


# --------------------------------------------------------------------------- #
# TSV
# --------------------------------------------------------------------------- #
def parse_tsv(raw: bytes | str) -> list[dict]:
    """解析样本 TSV，返回原始行列表。

    可选首行表头（首列不是整数即视为表头）。空行跳过。
    每行必须恰好 4 列：POS(1-based)、REF、ALT、GT。
    """
    text = _decode(raw)
    rows: list[dict] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        fields = [f.strip() for f in line.split("\t")]
        if not rows and fields and not _is_int(fields[0]):
            # 首行且首列不是数字 -> 视为表头
            continue
        if len(fields) != 4:
            raise ValidationError(f"第 {lineno} 行：应为 4 列（POS\\tREF\\tALT\\tGT），实际 {len(fields)} 列。")
        pos_s, ref, alt, gt = fields
        if not _is_int(pos_s):
            raise ValidationError(f"第 {lineno} 行：POS 必须是正整数，收到 '{pos_s}'。")
        rows.append(
            {"lineno": lineno, "pos": int(pos_s), "ref": ref, "alt": alt, "gt": gt}
        )
    if not rows:
        raise ValidationError("TSV 中没有任何变异记录。")
    return rows


def validate_and_normalize(rows: list[dict], ref_seq: str) -> list[dict]:
    """逐行核对 REF、校验类型并左对齐归一化。

    归一化后若同一样本出现重复变异（同一规范坐标），整份文件拒绝。
    返回 [{pos, ref, alt, gt, orig_pos, orig_ref, orig_alt, kind}, ...]。
    """
    genome_len = len(ref_seq)
    seen: dict[tuple, int] = {}
    result: list[dict] = []

    for row in rows:
        lineno = row["lineno"]
        pos, ref, alt, gt = row["pos"], row["ref"], row["alt"], row["gt"]

        if not ref or not alt:
            raise ValidationError(f"第 {lineno} 行：REF/ALT 不能为空。")
        if set(ref) - ACGT or set(alt) - ACGT:
            raise ValidationError(f"第 {lineno} 行：REF/ALT 只允许大写 A/C/G/T。")
        if gt not in VALID_GTS:
            raise ValidationError(f"第 {lineno} 行：GT 只允许 0/1 或 1/1，收到 '{gt}'。")
        if pos < 1:
            raise ValidationError(f"第 {lineno} 行：POS 必须从 1 开始计数。")

        kind = _classify(ref, alt, lineno)

        end = pos + len(ref) - 1
        if end > genome_len:
            raise ValidationError(
                f"第 {lineno} 行：变异结束位置 {end} 超出参考序列长度 {genome_len}。"
            )
        actual = ref_seq[pos - 1 : pos - 1 + len(ref)]
        if actual != ref:
            raise ValidationError(
                f"第 {lineno} 行：REF 核对失败，参考序列 {pos}-{end} 位为 '{actual}'，上传为 '{ref}'。"
            )

        if kind == "INS" and alt[0] != ref:
            raise ValidationError(f"第 {lineno} 行：插入的锚碱基必须与 REF 相同。")
        if kind == "DEL" and ref[0] != alt:
            raise ValidationError(f"第 {lineno} 行：缺失的锚碱基必须与 ALT 相同。")

        npos, nref, nalt = normalize_variant(pos, ref, alt, ref_seq)
        key = (npos, nref, nalt)
        if key in seen:
            raise ValidationError(
                f"第 {lineno} 行：归一化后与第 {seen[key]} 行为同一变异 "
                f"({npos}:{nref}>{nalt})，同一样本不允许重复变异，整份文件已拒绝。"
            )
        seen[key] = lineno

        result.append(
            {
                "pos": npos,
                "ref": nref,
                "alt": nalt,
                "gt": gt,
                "kind": kind,
                "orig_pos": pos,
                "orig_ref": ref,
                "orig_alt": alt,
            }
        )
    return result


def normalize_variant(pos: int, ref: str, alt: str, genome: str) -> tuple[int, str, str]:
    """把一个已通过校验的带锚变异左对齐到规范表示。

    左对齐过程中若在参考起点出现无法重新锚定的空等位基因
    （例如参考首位碱基上的纯重复插入），回退到最后一个合法的带锚状态，
    保证任何合法带锚输入都能得到合法的带锚规范表示。
    """
    # 单碱基替换不存在等价滑动，直接返回
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt

    last_anchored = (pos, ref, alt)
    while True:
        # 1) 去除共同后缀
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            if ref and alt:
                last_anchored = (pos, ref, alt)
            continue
        # 2) 无法去尾：左移并补参考前缀碱基
        if pos > 1:
            pos -= 1
            base = genome[pos - 1]
            ref = base + ref
            alt = base + alt
            last_anchored = (pos, ref, alt)
            continue
        # 3) 已到参考起点
        if ref and alt:
            break
        # 空等位基因且无法向左取锚：回退到最后一个合法带锚状态
        pos, ref, alt = last_anchored
        break

    # 4) 去除共同前缀，每个等位基因至少保留一个锚碱基
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1

    return pos, ref, alt


# --------------------------------------------------------------------------- #
# 内部工具
# --------------------------------------------------------------------------- #
def _classify(ref: str, alt: str, lineno: int) -> str:
    if len(ref) == 1 and len(alt) == 1:
        if ref == alt:
            raise ValidationError(f"第 {lineno} 行：REF 与 ALT 相同，不是变异。")
        return "SNP"
    if len(ref) == 1 and len(alt) >= 2:
        return "INS"
    if len(ref) >= 2 and len(alt) == 1:
        return "DEL"
    raise ValidationError(
        f"第 {lineno} 行：仅支持单碱基替换、带锚碱基的纯插入或纯缺失；"
        f"收到 REF='{ref}' ALT='{alt}'。"
    )


def _is_int(value: str) -> bool:
    return value.isdigit() and int(value) > 0


def _decode(raw: bytes | str) -> str:
    if isinstance(raw, bytes):
        try:
            return raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValidationError(f"文件不是合法的 UTF-8 文本：{exc}") from exc
    return raw
