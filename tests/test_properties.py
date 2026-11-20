"""随机化性质测试：用黄金性质（单倍型保持）验证归一化算法。

性质 1：任意合法变异归一化后，替换后单倍型与原始书写一致；
性质 2：两条变异若替换后单倍型相同，归一化键 (pos,ref,alt) 必须相同；
性质 3：在重复区内随机滑动重写的等价变异，归一化键必须相同。
"""

import random
import unittest

from app.normalization import normalize_variant

BASES = "ACGT"


def random_genome(rng, length):
    # 刻意制造低复杂度重复区
    if rng.random() < 0.5:
        unit = rng.choice(["A", "T", "AT", "TG", "GC"])
        return (unit * (length // len(unit) + 1))[:length]
    return "".join(rng.choice(BASES) for _ in range(length))


def haplotype(genome, pos, ref, alt):
    return genome[: pos - 1] + alt + genome[pos - 1 + len(ref):]


def random_variant(rng, genome):
    """在 genome 上生成一个合法的带锚 SNP/INS/DEL。"""
    pos = rng.randint(1, len(genome))
    kind = rng.choice(["SNP", "INS", "DEL"])
    if kind == "SNP":
        ref = genome[pos - 1]
        alt = rng.choice([b for b in BASES if b != ref])
        return pos, ref, alt
    if kind == "INS":
        ref = genome[pos - 1]
        n = rng.randint(1, 4)
        alt = ref + "".join(rng.choice(BASES) for _ in range(n))
        return pos, ref, alt
    # DEL：删除 pos 之后 1..3 个碱基（保证不越界）；末尾位点退化为 SNP
    n = rng.randint(1, min(3, len(genome) - pos)) if pos < len(genome) else 0
    if n == 0:
        ref = genome[pos - 1]
        alt = rng.choice([b for b in BASES if b != ref])
        return pos, ref, alt
    ref = genome[pos - 1 : pos + n]
    alt = genome[pos - 1]
    return pos, ref, alt


def slide_variant(rng, genome, pos, ref, alt):
    """尝试把变异在重复区内左右滑动，产生一个等价的不同书写；失败返回 None。"""
    for _ in range(10):
        delta = rng.choice([-1, 1])
        new_pos = pos + delta
        if new_pos < 1 or new_pos + len(ref) - 1 > len(genome):
            continue
        if delta < 1:
            base = genome[new_pos - 1]
            new_ref = base + ref[:-1]
            new_alt = base + alt[:-1] if len(alt) > 1 else base
        else:
            if new_pos + len(ref) - 1 > len(genome):
                continue
            base = genome[new_pos + len(ref) - 2]
            new_ref = ref[1:] + base
            new_alt = alt[1:] + base if len(alt) > 1 else alt
        # 必须仍然是参考支持的合法书写
        if genome[new_pos - 1 : new_pos - 1 + len(new_ref)] != new_ref:
            continue
        if len(new_ref) == 1 and len(new_alt) == 1 and new_ref == new_alt:
            continue
        return new_pos, new_ref, new_alt
    return None


class NormalizationPropertyTest(unittest.TestCase):
    def test_properties(self):
        rng = random.Random(20261119)
        # 键按「变异类型 + 长度变化量」分桶：同为缺失但删除碱基数不同的事件
        # 即使替换后单倍型碰巧相同，也不是同一个变异
        seen_haplotypes: dict[tuple, dict[tuple, tuple]] = {}
        checked_slides = 0

        for _ in range(4000):
            genome = random_genome(rng, rng.randint(8, 40))
            pos, ref, alt = random_variant(rng, genome)
            kind = "SNP" if len(ref) == 1 and len(alt) == 1 else (
                "INS" if len(ref) == 1 else "DEL"
            )
            event = (kind, abs(len(ref) - len(alt)))
            original_hap = haplotype(genome, pos, ref, alt)

            npos, nref, nalt = normalize_variant(pos, ref, alt, genome)
            # 性质 1：单倍型保持
            self.assertEqual(
                haplotype(genome, npos, nref, nalt),
                original_hap,
                f"haplotype changed: {(pos, ref, alt)} -> {(npos, nref, nalt)} on {genome}",
            )
            # 规范结果仍是合法带锚变异
            self.assertTrue(nref and nalt)
            self.assertGreaterEqual(npos, 1)
            # 变异类型与长度变化量在归一化前后保持不变
            nkind = "SNP" if len(nref) == 1 and len(nalt) == 1 else (
                "INS" if len(nref) == 1 else "DEL"
            )
            self.assertEqual((nkind, abs(len(nref) - len(nalt))), event)

            key = (npos, nref, nalt)
            # 性质 2：同一事件类型内，同单倍型 => 同规范键
            bucket = seen_haplotypes.setdefault(original_hap, {})
            if event in bucket:
                self.assertEqual(
                    bucket[event],
                    key,
                    f"same haplotype {original_hap} mapped to two {event} keys",
                )
            else:
                bucket[event] = key

            # 性质 3：等价滑动书写归一化到同一键
            slid = slide_variant(rng, genome, pos, ref, alt)
            if slid is not None:
                sp, sr, sa = slid
                if haplotype(genome, sp, sr, sa) == original_hap:
                    checked_slides += 1
                    self.assertEqual(
                        normalize_variant(sp, sr, sa, genome),
                        key,
                        f"equivalent slide diverged: {(pos, ref, alt)} vs {slid} on {genome}",
                    )

        self.assertGreater(checked_slides, 200, "随机滑动样本数不足，测试可信度不够")


if __name__ == "__main__":
    unittest.main()
