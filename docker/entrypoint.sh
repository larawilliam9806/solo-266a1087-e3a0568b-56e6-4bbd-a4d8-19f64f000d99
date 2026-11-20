#!/bin/sh
# 单 worker + 多线程：SQLite 写事务（BEGIN IMMEDIATE）天然串行化，
# 单进程避免跨进程锁竞争，线程内由 WAL + busy_timeout 保证一致性。
set -e

mkdir -p "$(dirname "$DATABASE_PATH")"

exec gunicorn --bind 0.0.0.0:8000 \
    --workers 1 \
    --threads 8 \
    --timeout 60 \
    --access-logfile - \
    --error-logfile - \
    "app:create_app()"
