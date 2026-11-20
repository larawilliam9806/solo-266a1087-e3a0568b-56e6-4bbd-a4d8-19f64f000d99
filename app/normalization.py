"""变异校验与参考左对齐归一化（纯函数，不依赖 Flask）。

支持的变异形式（均为 ACGT 大写、1-based 坐标）：
  * 单碱基替换：len(REF)==1 且 len(ALT)==1，REF != ALT
  * 带锚碱基的纯插入：len(REF)==1 且 len(ALT)>=2，且 ALT[0] == REF
  * 带锚碱基的纯缺失：len(REF)>=2 且 len(ALT)==1，且 REF[0] == ALT

归一化采用经典 left-align 流程（与 bcftools norm / vt normalize 等价）：
  1. 不断去掉两个等位基因相同的末尾碱基；
  2. 无法去尾时，把事件整体左移一位，在两个等位基因前补参考前缀碱基；
  3. 到达参考起点后，若出现空等位基因则回退到最后一个合法的带锚状态；
  4. 最后去除共同前缀（每条等位基因至少保留 1 个碱基）。
同一次插入/缺失在重复序列中无论以何种坐标书写，归一化结果都唯一。

性能：参考序列可能很长，而某条 indel 左侧大部分参考碱基与其事件无关。
逐位滑动的经典实现会让单条 indel 的耗时随无关前缀长度线性增长（字符串
反复拼接时甚至退化为平方复杂度）。这里用 :class:`BoundaryIndex` 保存参考
序列的连续相同碱基段（run）边界，把左对齐拆成两个阶段整体跳过：
  * 阶段一（步数不超过事件长度 n）：比较对象是变异自身的等位碱基，最多
    n 次，与参考长度完全无关；
  * 阶段二（步数超过 n）：退化为参考序列上「相距 n 的两个碱基是否相等」
    的周期判断，均一事件（n==1，最常见的单碱基插入/缺失）直接走 run
    边界 O(1) 定位，一般周期走指数倍增 + 二分，只访问与事件相关的区域。
因此单条 indel 的归一化计算量不随无关参考前缀的长度增长。
边界索引由调用方在「一次上传 / 一次换参考重验」开始时构建、结束即丢弃，
不做跨请求缓存（边界只在本次重验中有效）。
"""

from __future__ import annotations

from bisect import bisect_right

ACGT = frozenset("ACGT")
VALID_GTS = ("0/1", "1/1")
MAX_NAME_LEN = 128


class ValidationError(Exception):
    """输入数据不合法；message 可直接展示给用户。"""


# --------------------------------------------------------------------------- #
# 参考序列边界索引
# --------------------------------------------------------------------------- #
class BoundaryIndex:
    """参考序列的连续相同碱基段（run）边界索引。

    ``starts[i]`` 是第 i 个 run 的起始 0-based 位置，``bases[i]`` 是该 run
    的碱基。构建一次后可在同一次上传/重验的所有样本、所有变异间共用：
    多个样本若在同一段连续相同碱基中增减一个同类碱基，即使锚点不同，
    查的也是同一条 run 的同一条边界（结果按 (碱基, 位置) 记忆）。
    """

    __slots__ = ("genome", "starts", "bases", "_run_cache")

    def __init__(self, genome: str):
        self.genome = genome
        starts: list[int] = []
        bases: list[str] = []
        prev = ""
        for i, base in enumerate(genome):
            if base != prev:
                starts.append(i)
                bases.append(base)
                prev = base
        self.starts = starts
        self.bases = bases
        # (碱基, 0-based 位置) -> 所在 run 的起点；同 run 内只做一次二分
        self._run_cache: dict[tuple[str, int], int] = {}

    def run_start_at(self, pos0: int) -> int:
        """返回 0-based 位置 ``pos0`` 所在 run 的起始 0-based 位置。"""
        idx = bisect_right(self.starts, pos0) - 1
        return self.starts[idx]

    def homogeneous_run_start(self, base: str, pos0: int) -> int:
        """若 ``pos0`` 位于一段均为 ``base`` 的 run，返回该 run 的 0-based 起点。

        若该位置碱基不是 ``base``，返回 ``pos0``（长度为 0 的匹配）。
        """
        key = (base, pos0)
        cached = self._run_cache.get(key)
        if cached is not None:
            return cached
        if self.genome[pos0] != base:
            self._run_cache[key] = pos0
            return pos0
        start = self.run_start_at(pos0)
        self._run_cache[key] = start
        return start


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


def validate_and_normalize(
    rows: list[dict],
    ref_seq: str,
    index: BoundaryIndex | None = None,
) -> list[dict]:
    """逐行核对 REF、校验类型并左对齐归一化。

    归一化后若同一样本出现重复变异（同一规范坐标），整份文件拒绝。
    返回 [{pos, ref, alt, gt, orig_pos, orig_ref, orig_alt, kind}, ...]。

    ``index`` 为可选的 :class:`BoundaryIndex`，同一次换参考重验中由多个
    样本共用；不传则为本次调用单独构建一个（单样本上传路径）。
    """
    if index is None:
        index = BoundaryIndex(ref_seq)
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

        npos, nref, nalt = normalize_variant(pos, ref, alt, ref_seq, index)
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


def normalize_variant(
    pos: int,
    ref: str,
    alt: str,
    genome: str,
    index: BoundaryIndex | None = None,
) -> tuple[int, str, str]:
    """把一个已通过校验的带锚变异左对齐到规范表示。

    与逐位滑动的经典 left-align 完全等价（见模块 docstring），但通过
    :class:`BoundaryIndex` 按连续相同碱基段 / 周期块整体跳过：单条 indel
    的计算量只取决于与其事件相关的重复区域大小，不随左侧无关参考前缀
    的长度增长。
    """
    # 单碱基替换不存在等价滑动，直接返回
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt
    if index is None:
        index = BoundaryIndex(genome)

    # 1) 去除共同后缀
    stripped_ref, stripped_alt = ref, alt
    while stripped_ref and stripped_alt and stripped_ref[-1] == stripped_alt[-1]:
        stripped_ref = stripped_ref[:-1]
        stripped_alt = stripped_alt[:-1]

    # 去后缀后两条等位基因均非空：后续左移都会被结尾的共同前缀去除抵消，
    # 直接做共同前缀去除即可（与经典循环结果一致）
    if stripped_ref and stripped_alt:
        return _strip_common_prefix(pos, stripped_ref, stripped_alt)

    # 2) 去后缀后恰有一条等位基因变空：进入滑动事件
    p0 = pos
    max_slide = p0 - 1  # 最多左移到 1-based pos 1
    if max_slide == 0:
        # 已在参考起点且出现空等位基因：无法向左取锚，回退原带锚书写
        return _strip_common_prefix(p0, ref, alt)

    is_insertion = stripped_ref == ""
    # 去后缀后唯一非空的事件序列
    moving = stripped_alt if is_insertion else stripped_ref
    k = _left_slide_count(p0, moving, genome, index)

    # 3) 按滑动步数重建规范表示
    if k == max_slide:
        # 一直滑到参考起点仍为空等位基因：回退到最后一个合法带锚状态，
        # 即第 k 次左移之后、其共同后缀去除之前的状态。
        first_base = genome[0]
        tail = _rotate_right(moving, k - 1)
        if is_insertion:
            return _strip_common_prefix(1, first_base, first_base + tail)
        return _strip_common_prefix(1, first_base + tail, first_base)

    # 在参考起点之前遇到不匹配碱基而停下：第 k+1 次左移补入的碱基 B
    # （0-based 索引 p0-k-2），用它重新锚定；长等位为右移 k 次后的形态。
    left_base = genome[p0 - k - 2]
    moved = _rotate_right(moving, k)
    if is_insertion:
        return _strip_common_prefix(p0 - k - 1, left_base, left_base + moved)
    return _strip_common_prefix(p0 - k - 1, left_base + moved, left_base)


# --------------------------------------------------------------------------- #
# 内部工具
# --------------------------------------------------------------------------- #
def _left_slide_count(
    p0: int,
    moving: str,
    genome: str,
    index: BoundaryIndex,
) -> int:
    """计算空等位模式下可以连续成功左移的步数 k（k <= p0-1）。

    设去共同后缀后的非空事件序列为 ``moving``（长度 n）。第 t 次左移
    （t 从 1 开始）补入参考碱基 ``genome[p0-1-t]``，成功条件是该碱基
    等于事件序列右移 t-1 次后的末位，即：

        genome[p0-1-t] == moving[(n - t) mod n]

    分两个阶段求解，避免逐位扫描整段无关前缀：
      * 阶段一（t <= n）：比较对象是事件自身的碱基，最多 n 次；
      * 阶段二（t > n）：等式两边都是参考碱基且相距 n，即
        ``genome[i] == genome[i+n]``，按 run 边界 / 周期块整体跳过。
    """
    n = len(moving)
    limit = p0 - 1
    if limit == 0:
        return 0

    # 阶段一：t=1..min(n, limit)，期望序列（从右向左读）即 moving[::-1]
    period = moving[::-1]
    k1 = 0
    cap = min(n, limit)
    for t in range(1, cap + 1):
        if genome[p0 - 1 - t] == period[(t - 1) % n]:
            k1 += 1
        else:
            break
    if k1 >= limit:
        return limit
    if k1 < cap:
        return k1  # 阶段一内遇到不匹配，停止

    # 阶段二：t > n，条件 genome[i] == genome[i+n]
    return k1 + _periodic_equality_len(genome, index, n, p0, k1, limit)


def _periodic_equality_len(
    genome: str,
    index: BoundaryIndex,
    n: int,
    p0: int,
    k1: int,
    limit: int,
) -> int:
    """阶段二：统计连续满足 ``genome[i] == genome[i+n]`` 的步数。

    第一个待检查的左侧索引（t=k1+1 时）为 ``i0 = p0-1-(k1+1)``，随步数
    增加向左扩展，最多检查 ``limit-k1`` 个位置。
    """
    want = limit - k1
    if want <= 0:
        return 0
    i0 = p0 - 1 - (k1 + 1)

    if n == 1:
        # 单碱基事件：genome[i] == genome[i+1]，即 i 与其右邻同属一个 run。
        # 用 run 边界 O(1) 求出连续满足的个数，与无关前缀长度无关。
        if genome[i0] != genome[i0 + 1]:
            return 0
        run_start = index.homogeneous_run_start(genome[i0], i0)
        return min(want, i0 - run_start + 1)

    # 一般周期：指数倍增找上界，再二分求最长匹配；切片比较在 C 层完成，
    # 访问范围不超过真正匹配的区域 + 一块。
    def ok(m: int) -> bool:
        if m == 0:
            return True
        left = i0 - m + 1
        if left < 0:
            return False
        return genome[left : i0 + 1] == genome[left + n : i0 + 1 + n]

    hi = min(want, i0 + 1)
    while ok(hi) and hi < want:
        nxt = min(want, hi * 2)
        if nxt == hi:
            break
        hi = nxt
    lo = 0
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if ok(mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


def _rotate_right(seq: str, k: int) -> str:
    """序列循环右移 k 位。

    等价于经典流程中连续左移 k 次、每次在左侧补一个参考碱基后，
    事件序列的最终形态。
    """
    if not seq:
        return seq
    k %= len(seq)
    return seq[-k:] + seq[:-k] if k else seq


def _strip_common_prefix(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


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
