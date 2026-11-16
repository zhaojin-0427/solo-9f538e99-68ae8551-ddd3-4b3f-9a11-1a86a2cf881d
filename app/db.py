"""SQLite 访问层：连接管理、建表、中断任务恢复。"""
from __future__ import annotations

import os
import sqlite3

from flask import current_app, g

SCHEMA = """
CREATE TABLE IF NOT EXISTS batches (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    name        TEXT NOT NULL,
    created_at  TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS inputs (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id    INTEGER NOT NULL REFERENCES batches(id),
    kind        TEXT NOT NULL CHECK (kind IN ('fasta', 'csv')),
    version     INTEGER NOT NULL,
    content     TEXT NOT NULL,
    created_at  TEXT NOT NULL,
    UNIQUE (batch_id, kind, version)
);

CREATE TABLE IF NOT EXISTS jobs (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id        INTEGER NOT NULL REFERENCES batches(id),
    fasta_version   INTEGER NOT NULL,
    csv_version     INTEGER NOT NULL,
    status          TEXT NOT NULL DEFAULT 'pending'
                    CHECK (status IN ('pending', 'running', 'done', 'failed')),
    error           TEXT,
    created_at      TEXT NOT NULL,
    started_at      TEXT,
    finished_at     TEXT
);

-- 每个计算任务至多一份报告；报告整份写入（单事务），失败不留半份
CREATE TABLE IF NOT EXISTS reports (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id          INTEGER NOT NULL UNIQUE REFERENCES jobs(id),
    batch_id        INTEGER NOT NULL REFERENCES batches(id),
    fasta_version   INTEGER NOT NULL,
    csv_version     INTEGER NOT NULL,
    status          TEXT NOT NULL CHECK (status IN ('ok', 'failed')),
    error           TEXT,
    chrom           TEXT,
    ref_len         INTEGER,
    payload         TEXT,          -- ok 时为 JSON 结果；failed 时为 NULL
    created_at      TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status, id);
CREATE INDEX IF NOT EXISTS idx_reports_batch ON reports(batch_id, id);

-- 每个批次固定一份注释列定义（首个被接受的导入建立），之后列不一致整次拒绝
CREATE TABLE IF NOT EXISTS annotation_schemas (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id    INTEGER NOT NULL UNIQUE REFERENCES batches(id),
    columns     TEXT NOT NULL,          -- JSON 数组：注释列名（有序）
    created_at  TEXT NOT NULL
);

-- 注释按「批次 + 染色体名 + 完整参考序列哈希 + 标准位点」唯一存储：
-- 三者全部一致才可能被新报告复用，任一变化自然成为缺口，旧依据无法冒充新结论
CREATE TABLE IF NOT EXISTS annotations (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id        INTEGER NOT NULL REFERENCES batches(id),
    chrom           TEXT NOT NULL,
    ref_hash        TEXT NOT NULL,     -- 完整参考序列（大写）的 SHA256
    ref_len         INTEGER NOT NULL,
    pos             INTEGER NOT NULL,
    ref             TEXT NOT NULL,
    alt             TEXT NOT NULL,
    data            TEXT NOT NULL,     -- JSON：{注释列: 值}
    import_id       INTEGER NOT NULL REFERENCES annotation_imports(id),
    created_at      TEXT NOT NULL,
    UNIQUE (batch_id, chrom, ref_hash, pos, ref, alt)
);

CREATE INDEX IF NOT EXISTS idx_annotations_lookup
    ON annotations(batch_id, chrom, ref_hash);

-- 导入审计：成功（accepted）与失败（rejected）均留痕；
-- 失败导入不写入任何 annotations，不存在半次导入
CREATE TABLE IF NOT EXISTS annotation_imports (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id        INTEGER NOT NULL REFERENCES batches(id),
    report_id       INTEGER NOT NULL REFERENCES reports(id),
    filename        TEXT NOT NULL,
    status          TEXT NOT NULL CHECK (status IN ('accepted', 'rejected')),
    rows_seen       INTEGER NOT NULL DEFAULT 0,  -- 文件数据行数
    rows_inserted   INTEGER NOT NULL DEFAULT 0,  -- 实际新入库数（同值重传计为跳过）
    rows_skipped    INTEGER NOT NULL DEFAULT 0,  -- 同值重传跳过数
    error           TEXT,
    created_at      TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_annotation_imports_batch
    ON annotation_imports(batch_id, id);
"""


def connect(path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(path, timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=30000")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def get_db() -> sqlite3.Connection:
    if "db" not in g:
        g.db = connect(current_app.config["DATABASE"])
    return g.db


def close_db(_exc=None) -> None:
    db = g.pop("db", None)
    if db is not None:
        db.close()


def init_app(app) -> None:
    path = app.config["DATABASE"]
    if path != ":memory:":
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    conn = connect(path)
    try:
        conn.executescript(SCHEMA)
        # 进程重启后，此前 running 的任务视为被中断，重新排队
        conn.execute("UPDATE jobs SET status='pending', started_at=NULL WHERE status='running'")
        conn.commit()
    finally:
        conn.close()
    app.teardown_appcontext(close_db)
