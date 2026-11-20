"""归一化算法单元测试（标准库 unittest，无 Flask 依赖）。"""

import unittest

from app.normalization import (
    ValidationError,
    normalize_variant,
    parse_reference,
    parse_tsv,
    validate_and_normalize,
)


def canon(pos, ref, alt, genome):
    return normalize_variant(pos, ref, alt, genome)


class ParseReferenceTest(unittest.TestCase):
    def test_plain_text(self):
        self.assertEqual(parse_reference(b"ACGT\nACGT\n"), "ACGTACGT")

    def test_fasta_single_record(self):
        fasta = ">demo\nACGT\nTTTT\n"
        self.assertEqual(parse_reference(fasta), "ACGTTTTT")

    def test_fasta_multi_record_rejected(self):
        with self.assertRaises(ValidationError):
            parse_reference(">a\nACGT\n>b\nTTTT")

    def test_illegal_character(self):
        with self.assertRaises(ValidationError):
            parse_reference("ACGN")
        with self.assertRaises(ValidationError):
            parse_reference("acgt")  # 小写不自动转换

    def test_empty(self):
        with self.assertRaises(ValidationError):
            parse_reference("")
        with self.assertRaises(ValidationError):
            parse_reference(">only header\n")


class ParseTsvTest(unittest.TestCase):
    def test_with_header(self):
        rows = parse_tsv("POS\tREF\tALT\tGT\n3\tA\tT\t0/1\n")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["pos"], 3)
        self.assertEqual(rows[0]["gt"], "0/1")
        self.assertEqual(rows[0]["lineno"], 2)

    def test_without_header(self):
        rows = parse_tsv("3\tA\tT\t0/1\n5\tC\tG\t1/1\n")
        self.assertEqual(len(rows), 2)

    def test_blank_lines_skipped(self):
        rows = parse_tsv("\n3\tA\tT\t0/1\n\n")
        self.assertEqual(len(rows), 1)

    def test_wrong_column_count(self):
        with self.assertRaises(ValidationError):
            parse_tsv("3\tA\tT\n")

    def test_bad_pos(self):
        with self.assertRaises(ValidationError):
            parse_tsv("x\tA\tT\t0/1\n")
        with self.assertRaises(ValidationError):
            parse_tsv("0\tA\tT\t0/1\n")

    def test_empty_file(self):
        with self.assertRaises(ValidationError):
            parse_tsv("POS\tREF\tALT\tGT\n")


class SnpTest(unittest.TestCase):
    GENOME = "ACGTACGT"

    def test_snp_unchanged(self):
        self.assertEqual(canon(3, "G", "T", self.GENOME), (3, "G", "T"))

    def test_same_ref_alt_rejected(self):
        with self.assertRaises(ValidationError):
            validate_and_normalize(
                [{"lineno": 1, "pos": 1, "ref": "A", "alt": "A", "gt": "0/1"}],
                self.GENOME,
            )


class IndelEquivalenceTest(unittest.TestCase):
    """核心：重复序列中不同坐标书写的同一 indel 必须归一化为同一规范变异。"""

    # 位置:    1 2 3 4 5 6 7
    GENOME = "G C T T T G C".replace(" ", "")  # GCTTTGC

    def _haplotype(self, pos, ref, alt):
        """变异对应的替换后单倍型（用于等价性验证，与书写坐标无关）。"""
        return self.GENOME[: pos - 1] + alt + self.GENOME[pos - 1 + len(ref):]

    def test_insertion_equivalent_forms(self):
        # 在 TTT 重复区插入一个 T，三个锚定位点书写等价；
        # 最左可滑动到 C 锚（C 后接 T 与插入碱基相同），规范为 (2,C,CT)
        forms = [(3, "T", "TT"), (4, "T", "TT"), (5, "T", "TT")]
        results = {canon(p, r, a, self.GENOME) for p, r, a in forms}
        self.assertEqual(results, {(2, "C", "CT")})
        # 三种书写的替换后单倍型必须一致
        haps = {self._haplotype(p, r, a) for p, r, a in forms}
        self.assertEqual(len(haps), 1)

    def test_insertion_non_equivalent_forms_differ(self):
        # 在 pos1 G 后插入 C，与在重复区插入 T 不是同一变异
        self.assertNotEqual(
            canon(1, "G", "GC", self.GENOME),
            canon(3, "T", "TT", self.GENOME),
        )

    def test_deletion_equivalent_forms(self):
        # 删除 TTT 中同一个 T：(3,TT,T) 与 (4,TT,T) 等价，规范为 (2,CT,C)
        forms = [(3, "TT", "T"), (4, "TT", "T")]
        results = {canon(p, r, a, self.GENOME) for p, r, a in forms}
        self.assertEqual(results, {(2, "CT", "C")})
        haps = {self._haplotype(p, r, a) for p, r, a in forms}
        self.assertEqual(haps, {"GCTTGC"})  # GCTTTGC 删掉一个 T

    def test_different_haplotypes_not_merged(self):
        # 单倍型不同的缺失写法不能被错误合并
        # (5,TT,T): GCTTTC；(3,GTT,GT): GCGTGC
        self.assertNotEqual(
            canon(5, "TT", "T", self.GENOME),
            canon(3, "GTT", "GT", self.GENOME),
        )
        self.assertNotEqual(
            self._haplotype(5, "TT", "T"),
            self._haplotype(3, "GTT", "GT"),
        )

    def test_longer_anchor_forms_equivalent(self):
        # 同一删除（单倍型 GCTTGC，删 pos3 的 T）的短锚/长锚/规范书写三者等价
        forms = [(3, "TT", "T"), (4, "TT", "T"), (2, "CT", "C")]
        results = {canon(p, r, a, self.GENOME) for p, r, a in forms}
        self.assertEqual(len(results), 1)
        self.assertEqual(results.pop(), (2, "CT", "C"))

    def test_different_deletion_lengths_are_not_equivalent(self):
        # 不同单倍型的缺失不能因为锚的重叠而被错误合并
        distinct = [
            (3, "GTT", "GT"),   # GCGTGC
            (5, "TT", "T"),     # GCTTTC
            (3, "GCTT", "GCT"), # GCGCTC
        ]
        canon_set = {canon(p, r, a, self.GENOME) for p, r, a in distinct}
        self.assertEqual(len(canon_set), len(distinct))

    def test_haplotype_preserved_insertion(self):
        # 归一化前后单倍型含义必须保持一致
        original_hap = self._haplotype(5, "T", "TT")
        p, r, a = canon(5, "T", "TT", self.GENOME)
        self.assertEqual(self._haplotype(p, r, a), original_hap)

    def test_haplotype_preserved_deletion(self):
        original_hap = self._haplotype(3, "TT", "T")
        p, r, a = canon(3, "TT", "T", self.GENOME)
        self.assertEqual(self._haplotype(p, r, a), original_hap)

    def test_leftmost_position_is_canonical(self):
        p, _, _ = canon(5, "T", "TT", self.GENOME)
        self.assertEqual(p, 2)
        p, _, _ = canon(5, "TT", "T", self.GENOME)
        self.assertEqual(p, 2)


class BoundaryTest(unittest.TestCase):
    def test_insertion_at_position_one(self):
        genome = "ACGT"
        self.assertEqual(canon(1, "A", "AT", genome), (1, "A", "AT"))

    def test_deletion_at_position_one(self):
        # CAAA 上 (2,AA,A) 左对齐到 (1,CA,C)
        genome = "CAAA"
        self.assertEqual(canon(2, "AA", "A", genome), (1, "CA", "C"))

    def test_insertion_reanchored_at_left_edge(self):
        # pos=1 插入与首碱基不同的碱基：无法左移，保持原位
        genome = "ACGT"
        self.assertEqual(canon(1, "A", "AC", genome), (1, "A", "AC"))

    def test_insertion_reanchored_via_shift(self):
        # CAACG 上 pos=2 的 A 插入向左滑到 pos=1，空等位基因向左取 C 锚重新锚定
        genome = "CAACG"
        self.assertEqual(canon(2, "A", "AA", genome), (1, "C", "CA"))
        # pos=1 的等价书写归一化结果相同，且单倍型一致
        self.assertEqual(canon(1, "C", "CA", genome), (1, "C", "CA"))
        p, ref, alt = canon(2, "A", "AA", genome)
        self.assertEqual(
            genome[: p - 1] + alt + genome[p - 1 + len(ref):],
            "CAAACG",
        )

    def test_unanchorable_at_position_one_keeps_last_anchored(self):
        # pos=1 插入与首碱基相同的碱基，左移到起点后无法向左取锚，
        # 回退到最后一个合法带锚状态（即原书写本身）
        genome = "ACGT"
        self.assertEqual(canon(1, "A", "AA", genome), (1, "A", "AA"))

    def test_unanchorable_deletion_keeps_last_anchored(self):
        # AAAC 上删除首个 A 的等价形式滑到起点后无法锚定，回退到合法带锚状态
        genome = "AAAC"
        result = canon(2, "AA", "A", genome)
        p, ref, alt = result
        self.assertTrue(ref and alt)
        self.assertEqual(
            genome[: p - 1] + alt + genome[p - 1 + len(ref):],
            "AAC",
        )

    def test_out_of_range_rejected(self):
        genome = "ACGT"
        with self.assertRaises(ValidationError):
            validate_and_normalize(
                [{"lineno": 1, "pos": 4, "ref": "TT", "alt": "T", "gt": "0/1"}],
                genome,
            )


class BatchValidationTest(unittest.TestCase):
    GENOME = "GCTTTGC"

    def _row(self, pos, ref, alt, gt="0/1", lineno=1):
        return {"lineno": lineno, "pos": pos, "ref": ref, "alt": alt, "gt": gt}

    def test_ref_mismatch_rejected(self):
        with self.assertRaises(ValidationError) as ctx:
            validate_and_normalize([self._row(1, "C", "G")], self.GENOME)
        self.assertIn("REF 核对失败", str(ctx.exception))

    def test_unsupported_variant_type(self):
        with self.assertRaises(ValidationError):
            # 替换 + 插入混合（REF/ALT 都长且不等长），不支持
            validate_and_normalize([self._row(1, "AC", "GT")], self.GENOME)

    def test_bad_gt(self):
        with self.assertRaises(ValidationError):
            validate_and_normalize([self._row(1, "G", "C", gt="0/0")], self.GENOME)
        with self.assertRaises(ValidationError):
            validate_and_normalize([self._row(1, "G", "C", gt="1|1")], self.GENOME)

    def test_anchor_mismatch_rejected(self):
        with self.assertRaises(ValidationError):
            # 插入锚碱基必须等于 REF
            validate_and_normalize([self._row(1, "G", "CT")], self.GENOME)
        with self.assertRaises(ValidationError):
            # 缺失锚碱基必须等于 ALT
            validate_and_normalize([self._row(1, "GC", "T")], self.GENOME)

    def test_normalized_duplicate_rejects_whole_file(self):
        rows = [
            self._row(3, "T", "TT", lineno=1),
            self._row(5, "T", "TT", lineno=2),  # 归一化后与第一行重复
        ]
        with self.assertRaises(ValidationError) as ctx:
            validate_and_normalize(rows, self.GENOME)
        self.assertIn("同一变异", str(ctx.exception))

    def test_distinct_variants_pass(self):
        rows = [
            self._row(1, "G", "C", lineno=1, gt="1/1"),
            self._row(3, "T", "TT", lineno=2, gt="0/1"),
        ]
        out = validate_and_normalize(rows, self.GENOME)
        self.assertEqual(len(out), 2)
        self.assertEqual(out[0]["kind"], "SNP")
        self.assertEqual(out[1]["kind"], "INS")
        self.assertEqual(out[1]["orig_pos"], 3)

    def test_exact_duplicate_rows_rejected(self):
        rows = [
            self._row(3, "T", "TT", lineno=1),
            self._row(3, "T", "TT", lineno=2),
        ]
        with self.assertRaises(ValidationError):
            validate_and_normalize(rows, self.GENOME)


if __name__ == "__main__":
    unittest.main()
