"""端到端集成测试：Flask 测试客户端 + 临时 SQLite。"""

import io
import os
import tempfile
import unittest

from app import create_app, services

REF_V1 = "GCTTTGCACGT"          # 11 bp，含 TTT 重复区
REF_V2 = "GCTATGCACGT"          # 第 4 位 T->A，原 TTT 区结构变化


def tsv(text: str) -> bytes:
    return text.encode("utf-8")


class AppTestBase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db_path = os.path.join(self.tmp.name, "test.db")
        os.environ["DATABASE_PATH"] = self.db_path
        self.app = create_app()
        self.app.config.update(TESTING=True, SECRET_KEY="test")
        self.client = self.app.test_client()

    def tearDown(self):
        os.environ.pop("DATABASE_PATH", None)
        self.tmp.cleanup()

    # --- 辅助 ------------------------------------------------------------- #
    def upload_ref(self, seq: bytes | str):
        payload = seq.encode() if isinstance(seq, str) else seq
        data = {"reference": (io.BytesIO(payload), "ref.txt")}
        return self.client.post(
            "/reference", data=data, content_type="multipart/form-data",
            follow_redirects=True,
        )

    def upload_sample(self, name: str, content: bytes):
        data = {"name": name, "sample": (io.BytesIO(content), f"{name}.tsv")}
        return self.client.post(
            "/samples", data=data, content_type="multipart/form-data",
            follow_redirects=True,
        )

    def db_query(self, sql, args=()):
        with self.app.app_context():
            from app.db import get_db
            return get_db().execute(sql, args).fetchall()

    def dashboard_map(self):
        with self.app.app_context():
            dash = services.get_dashboard()
        return {
            (v["pos"], v["ref"], v["alt"]): v["n_samples"]
            for v in dash["variants"]
        }


class ReferenceFirstTest(AppTestBase):
    def test_sample_before_reference_rejected(self):
        resp = self.upload_sample("S1", tsv("1\tG\tC\t0/1\n"))
        self.assertIn("请先上传参考序列", resp.get_data(as_text=True))
        self.assertEqual(self.db_query("SELECT COUNT(*) c FROM samples")[0]["c"], 0)

    def test_bad_reference_rejected(self):
        resp = self.upload_ref("ACGN")
        self.assertIn("只允许包含大写", resp.get_data(as_text=True))
        self.assertEqual(self.db_query("SELECT COUNT(*) c FROM reference_sequences")[0]["c"], 0)

    def test_fasta_accepted(self):
        resp = self.upload_ref(">demo ref\nACGT\nACGT\n")
        self.assertIn("参考已更新", resp.get_data(as_text=True))
        row = self.db_query("SELECT length FROM reference_sequences")[0]
        self.assertEqual(row["length"], 8)


class FullWorkflowTest(AppTestBase):
    def setUp(self):
        super().setUp()
        self.upload_ref(REF_V1)

    def test_equivalent_variants_merged_across_samples(self):
        # S1: 重复区内 pos=3 书写插入 T；另含一个 SNP
        self.upload_sample(
            "S1", tsv("POS\tREF\tALT\tGT\n3\tT\tTT\t0/1\n8\tA\tC\t1/1\n")
        )
        # S2: 同一插入用 pos=5 书写（等价），SNP 同一位点、不同 GT
        self.upload_sample("S2", tsv("5\tT\tTT\t0/1\n8\tA\tC\t0/1\n"))
        # S3: 等价缺失
        self.upload_sample("S3", tsv("3\tTT\tT\t1/1\n"))

        summary = self.dashboard_map()
        # 两种坐标的插入归一化为同一规范变异 (2,C,CT)，2 个样本检出
        self.assertEqual(summary[(2, "C", "CT")], 2)
        # SNP 2 个样本
        self.assertEqual(summary[(8, "A", "C")], 2)
        # 缺失规范为 (2,CT,C)，1 个样本
        self.assertEqual(summary[(2, "CT", "C")], 1)
        self.assertEqual(len(summary), 3)

        # 详情页保留各自原始位点与 GT
        detail = self.client.get("/variants/2/C>CT").get_data(as_text=True)
        self.assertIn("S1", detail)
        self.assertIn("S2", detail)
        self.assertIn(">3<", detail)   # S1 原始 POS=3
        self.assertIn(">5<", detail)   # S2 原始 POS=5
        self.assertNotIn("S3", detail)

    def test_duplicate_after_normalization_rejects_whole_sample(self):
        resp = self.upload_sample("BAD", tsv("3\tT\tTT\t0/1\n5\tT\tTT\t0/1\n"))
        self.assertIn("归一化后与第 1 行为同一变异", resp.get_data(as_text=True))
        self.assertEqual(self.db_query("SELECT COUNT(*) c FROM samples")[0]["c"], 0)
        self.assertEqual(
            self.db_query("SELECT COUNT(*) c FROM sample_variants")[0]["c"], 0
        )

    def test_ref_mismatch_rejected_and_state_unchanged(self):
        self.upload_sample("S1", tsv("1\tG\tC\t0/1\n"))
        resp = self.upload_sample("BAD", tsv("1\tT\tC\t0/1\n"))
        self.assertIn("REF 核对失败", resp.get_data(as_text=True))
        # 旧样本与旧变异不受影响
        self.assertEqual(self.db_query("SELECT COUNT(*) c FROM samples")[0]["c"], 1)
        self.assertEqual(
            self.db_query("SELECT COUNT(*) c FROM sample_variants")[0]["c"], 1
        )

    def test_invalid_gt_rejected(self):
        for bad_gt in ("0/0", "1|1", "./.", "1/0"):
            resp = self.upload_sample("BAD", tsv(f"1\tG\tC\t{bad_gt}\n"))
            self.assertIn("GT 只允许", resp.get_data(as_text=True))
        self.assertEqual(self.db_query("SELECT COUNT(*) c FROM samples")[0]["c"], 0)

    def test_unsupported_variant_type_rejected(self):
        # 多碱基替换（REF/ALT 均长且不等长规则之外）不接收
        resp = self.upload_sample("BAD", tsv("1\tGC\tAT\t0/1\n"))
        self.assertIn("仅支持", resp.get_data(as_text=True))
        # 锚碱基不匹配的缺失不接收（REF 与参考一致但首碱基 != ALT）
        resp = self.upload_sample("BAD", tsv("6\tGC\tC\t0/1\n"))
        self.assertIn("锚碱基", resp.get_data(as_text=True))
        # 多碱基的带锚纯插入是合法的，应通过
        resp = self.upload_sample("OKINS", tsv("1\tG\tGCAT\t0/1\n"))
        self.assertIn("成功", resp.get_data(as_text=True))

    def test_reupload_replaces_old_version(self):
        self.upload_sample("S1", tsv("1\tG\tC\t0/1\n8\tA\tC\t1/1\n"))
        self.assertEqual(
            self.db_query("SELECT COUNT(*) c FROM sample_variants")[0]["c"], 2
        )
        resp = self.upload_sample("S1", tsv("1\tG\tC\t1/1\n"))
        self.assertIn("替换旧版本成功", resp.get_data(as_text=True))
        rows = self.db_query("SELECT * FROM sample_variants")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["gt"], "1/1")
        self.assertEqual(
            self.db_query("SELECT COUNT(*) c FROM samples WHERE name='S1'")[0]["c"], 1
        )

    def test_failed_reupload_keeps_old_version(self):
        self.upload_sample("S1", tsv("1\tG\tC\t0/1\n"))
        resp = self.upload_sample("S1", tsv("1\tT\tC\t0/1\n"))  # REF 不符
        self.assertIn("上传失败", resp.get_data(as_text=True))
        rows = self.db_query(
            "SELECT s.name, s.status, v.gt FROM samples s "
            "JOIN sample_variants v ON v.sample_id=s.id"
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["gt"], "0/1")
        self.assertEqual(rows[0]["status"], "VALID")

    def test_reference_change_revalidation(self):
        self.upload_sample("S1", tsv("1\tG\tC\t0/1\n3\tT\tTT\t0/1\n"))
        # S2 依赖 V1 第 4 位的 T
        self.upload_sample("S2", tsv("4\tT\tTT\t1/1\n"))

        resp = self.upload_ref(REF_V2)
        self.assertIn("参考已更新", resp.get_data(as_text=True))

        statuses = {
            r["name"]: r["status"]
            for r in self.db_query("SELECT name, status FROM samples")
        }
        # S1 的变异在 V2 下 REF 仍匹配，重验通过
        self.assertEqual(statuses["S1"], "VALID")
        # S2 第 4 位在 V2 下为 A，REF 不符 -> 待重传并从汇总剔除
        self.assertEqual(statuses["S2"], "NEEDS_REUPLOAD")

        counts = {
            r["name"]: r["c"]
            for r in self.db_query(
                "SELECT s.name, COUNT(v.id) c FROM samples s "
                "LEFT JOIN sample_variants v ON v.sample_id=s.id GROUP BY s.name"
            )
        }
        self.assertEqual(counts["S2"], 0)
        self.assertGreater(counts["S1"], 0)

        page = self.client.get("/").get_data(as_text=True)
        self.assertIn("待重传", page)
        self.assertIn("REF 核对失败", page)
        # 汇总中只剩 S1 贡献的变异
        for n in self.dashboard_map().values():
            self.assertEqual(n, 1)

        # S2 按新参考修正重传后恢复
        resp = self.upload_sample("S2", tsv("4\tA\tAT\t1/1\n"))
        self.assertIn("成功", resp.get_data(as_text=True))
        row = self.db_query("SELECT status FROM samples WHERE name='S2'")[0]
        self.assertEqual(row["status"], "VALID")

    def test_summary_only_counts_current_reference(self):
        self.upload_sample("S1", tsv("1\tG\tC\t0/1\n"))
        self.upload_sample("S2", tsv("1\tG\tC\t1/1\n"))
        self.assertEqual(self.dashboard_map()[(1, "G", "C")], 2)
        # 换参考后，两个样本都被重验并指向新参考；S1/S2 的 SNP 在 V2 第 1 位仍为 G
        self.upload_ref(REF_V2)
        self.assertEqual(self.dashboard_map()[(1, "G", "C")], 2)
        page = self.client.get("/").get_data(as_text=True)
        self.assertNotIn("历史参考", page)

    def test_delete_sample_updates_counts(self):
        self.upload_sample("S1", tsv("1\tG\tC\t0/1\n"))
        self.upload_sample("S2", tsv("1\tG\tC\t1/1\n"))
        sid = self.db_query("SELECT id FROM samples WHERE name='S1'")[0]["id"]
        self.client.post(f"/samples/{sid}/delete", follow_redirects=True)
        self.assertEqual(
            self.db_query("SELECT COUNT(*) c FROM samples")[0]["c"], 1
        )
        detail = self.client.get("/variants/1/G>C").get_data(as_text=True)
        self.assertIn("S2", detail)
        self.assertNotIn("S1", detail)

    def test_healthz_and_404(self):
        self.assertEqual(self.client.get("/healthz").get_json(), {"status": "ok"})
        self.assertEqual(self.client.get("/variants/999/A>T").status_code, 404)

    def test_missing_file_fields(self):
        resp = self.client.post(
            "/samples", data={"name": "x"}, content_type="multipart/form-data",
            follow_redirects=True,
        )
        self.assertIn("请选择样本 TSV 文件", resp.get_data(as_text=True))


if __name__ == "__main__":
    unittest.main()
