"""SQLite 连接管理与事务。

并发模型：
  * 开启 WAL，允许多读单写；写事务一律 BEGIN IMMEDIATE 立即取写锁，
    使「样本上传」与「参考更新」严格串行化，网页读到的永远是
    某一版参考下完整、一致的汇总。
  * Python sqlite3 默认会在 DML 前自动 BEGIN，这里显式关闭
    （isolation_level=None），事务边界完全由本模块控制。
"""

from __future__ import annotations

import sqlite3
from contextlib import contextmanager

from flask import Flask, current_app, g


def get_db() -> sqlite3.Connection:
    if "db" not in g:
        conn = sqlite3.connect(
            current_app.config["DATABASE"],
            detect_types=sqlite3.PARSE_DECLTYPES,
        )
        conn.row_factory = sqlite3.Row
        conn.isolation_level = None  # 手动事务
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("PRAGMA busy_timeout=5000")
        g.db = conn
    return g.db


def close_db(_exc=None) -> None:
    conn: sqlite3.Connection | None = g.pop("db", None)
    if conn is not None:
        conn.close()


@contextmanager
def read_tx():
    """读事务：单条快照，保证一次请求内所有读取来自同一版本。"""
    conn = get_db()
    conn.execute("BEGIN")
    try:
        yield conn
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise


@contextmanager
def write_tx():
    """写事务：立即取写锁，串行化所有上传/更新操作；失败整体回滚。"""
    conn = get_db()
    conn.execute("BEGIN IMMEDIATE")
    try:
        yield conn
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise


SCHEMA = """
CREATE TABLE IF NOT EXISTS reference_sequences (
    id            INTEGER PRIMARY KEY,
    sequence      TEXT NOT NULL,
    length        INTEGER NOT NULL,
    uploaded_at   TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS samples (
    id                INTEGER PRIMARY KEY,
    name              TEXT NOT NULL UNIQUE,
    reference_id      INTEGER REFERENCES reference_sequences(id),
    status            TEXT NOT NULL CHECK (status IN ('VALID', 'NEEDS_REUPLOAD')),
    error             TEXT,
    raw_tsv           TEXT NOT NULL,
    uploaded_at       TEXT NOT NULL DEFAULT (datetime('now')),
    revalidated_at    TEXT
);

CREATE TABLE IF NOT EXISTS sample_variants (
    id            INTEGER PRIMARY KEY,
    sample_id     INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
    pos           INTEGER NOT NULL,
    ref           TEXT NOT NULL,
    alt           TEXT NOT NULL,
    gt            TEXT NOT NULL,
    kind          TEXT NOT NULL,
    orig_pos      INTEGER NOT NULL,
    orig_ref      TEXT NOT NULL,
    orig_alt      TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_variant_canonical
    ON sample_variants(pos, ref, alt);
CREATE INDEX IF NOT EXISTS idx_variant_sample ON sample_variants(sample_id);
"""


def init_db(app: Flask) -> None:
    conn = sqlite3.connect(app.config["DATABASE"])
    try:
        conn.executescript(SCHEMA)
        conn.commit()
    finally:
        conn.close()


def init_app(app: Flask) -> None:
    app.teardown_appcontext(close_db)
    with app.app_context():
        init_db(app)
