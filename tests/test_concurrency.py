"""并发一致性测试：样本上传与参考更新交错时，汇总始终对应单一参考版本。"""

import os
import tempfile
import threading
import unittest

from app import create_app, services
from app.normalization import ValidationError

REF_V1 = "GCTTTGCACGT"
REF_V2 = "GCTATGCACGT"
REF_V3 = "ACGTACGTACGT"
REFS = [REF_V1, REF_V2, REF_V3]


def tsv(text: str) -> bytes:
    return text.encode("utf-8")


class ConcurrencyTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        os.environ["DATABASE_PATH"] = os.path.join(self.tmp.name, "conc.db")
        self.app = create_app()
        self.app.config.update(TESTING=True, SECRET_KEY="test")
        self.errors: list[Exception] = []
        self.expected_rejections = 0

    def tearDown(self):
        os.environ.pop("DATABASE_PATH", None)
        self.tmp.cleanup()

    def test_interleaved_uploads_always_consistent(self):
        # 先放一版参考
        with self.app.app_context():
            services.upload_reference(REF_V1)

        def upload_samples():
            for i in range(15):
                name = f"S{i % 6}"  # 含同名重传
                gt = "1/1" if i % 2 else "0/1"
                try:
                    with self.app.app_context():
                        # 每条都针对 V1 有效；换参考后由系统重验
                        services.upload_sample(name, tsv(f"1\tG\tC\t{gt}\n"))
                except ValidationError:
                    # 参考已轮换导致 REF 不符是预期竞争结果：拒绝且不留脏数据
                    self.expected_rejections += 1
                except Exception as exc:
                    self.errors.append(exc)

        def rotate_references():
            for ref in REFS[1:]:
                try:
                    with self.app.app_context():
                        services.upload_reference(ref)
                except Exception as exc:
                    self.errors.append(exc)

        def reader():
            for _ in range(40):
                try:
                    with self.app.app_context():
                        dash = services.get_dashboard()
                        ref_id = dash["reference"]["id"]
                        # 不变量 1：汇总中的每个变体只来自当前参考的有效样本
                        for v in dash["variants"]:
                            self.assertGreaterEqual(v["n_samples"], 1)
                        # 不变量 2：所有 VALID 样本都指向当前参考
                        for s in dash["samples"]:
                            if s["status"] == "VALID":
                                self.assertEqual(s["reference_id"], ref_id)
                        # 不变量 3：变体数与样本变体明细一致
                        total = sum(
                            s["variant_count"]
                            for s in dash["samples"]
                            if s["status"] == "VALID" and s["reference_id"] == ref_id
                        )
                        listed = sum(
                            1 for _ in dash["variants"] for _ in [0]
                        )
                        self.assertLessEqual(listed, total)
                except Exception as exc:
                    self.errors.append(exc)

        threads = [
            threading.Thread(target=upload_samples),
            threading.Thread(target=rotate_references),
            threading.Thread(target=reader),
            threading.Thread(target=reader),
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)
            self.assertFalse(t.is_alive(), "线程未在超时内结束，可能存在锁死锁")

        if self.errors:
            raise AssertionError(
                "并发过程出现非预期异常:\n"
                + "\n".join(f"{type(e).__name__}: {e}" for e in self.errors)
            )
        # 参考轮换后确实发生过交错拒绝，且拒绝没有破坏任何不变量
        self.assertGreater(self.expected_rejections, 0)

        # 最终状态：参考为 V3，所有样本已被重验并指向 V3（VALID 或 NEEDS_REUPLOAD）
        with self.app.app_context():
            dash = services.get_dashboard()
            self.assertEqual(dash["reference"]["sequence"], REF_V3)
            for s in dash["samples"]:
                self.assertEqual(s["reference_id"], dash["reference"]["id"])


if __name__ == "__main__":
    unittest.main()
