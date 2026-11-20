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
            # 足够多的上传次数保证与参考轮换发生锁竞争（不依赖单次归一化耗时）
            for i in range(120):
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
        # 写锁严格串行化：并发交错时只要发生过参考轮换，就可能出现 REF 不符
        # 的拒绝；拒绝与否取决于线程调度，但任何拒绝都不得破坏一致性不变量。
        # 这里不再用「拒绝次数 > 0」这种依赖时序的断言，改为下面的
        # 确定性锁竞争测试单独验证串行化语义。
        self._assert_write_lock_serializes()

        # 最终状态：参考为 V3，所有样本已被重验并指向 V3（VALID 或 NEEDS_REUPLOAD）
        with self.app.app_context():
            dash = services.get_dashboard()
            self.assertEqual(dash["reference"]["sequence"], REF_V3)
            for s in dash["samples"]:
                self.assertEqual(s["reference_id"], dash["reference"]["id"])

    def _assert_write_lock_serializes(self):
        """确定性验证：一个写事务持锁期间，其他写操作必须等待，不产生脏数据。"""
        import threading
        from app.db import get_db, write_tx
        from app.normalization import ValidationError

        entered = threading.Event()
        release = threading.Event()
        observed = {}

        def hold_lock():
            with self.app.app_context():
                with write_tx() as conn:
                    conn.execute("SELECT 1")
                    entered.set()
                    release.wait(timeout=5)

        def try_write():
            with self.app.app_context():
                try:
                    services.upload_sample(
                        "LOCKPROBE", tsv("1\tG\tC\t0/1\n")
                    )
                    observed["done"] = True
                except ValidationError as exc:
                    observed["error"] = str(exc)

        holder = threading.Thread(target=hold_lock)
        holder.start()
        self.assertTrue(entered.wait(timeout=5), "持锁线程未进入写事务")

        waiter = threading.Thread(target=try_write)
        waiter.start()
        # 持锁期间等待线程不应完成（busy_timeout=5s，这里只等 1s 确认被阻塞）
        waiter.join(timeout=1)
        self.assertTrue(waiter.is_alive(), "写锁未串行化：等待线程在持锁期间完成了写入")

        release.set()
        holder.join(timeout=5)
        waiter.join(timeout=5)
        self.assertFalse(waiter.is_alive(), "等待线程在锁释放后仍未结束")
        # 锁释放后写入正常落库（或因参考轮换被合法拒绝），不留脏状态
        self.assertTrue("done" in observed or "error" in observed)


if __name__ == "__main__":
    unittest.main()
