"""注释 CSV 导入、覆盖情况与按源 ID 展开下载。

复用语义（均由查询条件派生，不做任何状态标记）：
- 一条注释只对「同批次 + 同染色体名 + 同完整参考序列哈希 + 同标准位点」可见；
  参考更正或位点变化后，旧注释自动成为新报告的缺口，绝不可能冒充当前完整结果；
- 导入与输入更正/后台重算以任意顺序交错都安全：可复用注释按依据即时计算，
  既不依赖任务先后，也不会漏掉满足依据的历史注释。

原子性：一次导入在单事务（BEGIN IMMEDIATE）内完成，任何错误整体回滚，
仅留一条 rejected 审计记录，不留半次导入。
"""
from __future__ import annotations

import json

from flask import (Blueprint, Response, abort, flash, redirect,
                   render_template, request, url_for)

from . import annotate
from .db import get_db
from .normalize import ValidationError

bp = Blueprint("annotations", __name__)


# ---------------------------------------------------------------- 辅助

def _batch_or_404(db, bid):
    batch = db.execute("SELECT * FROM batches WHERE id=?", (bid,)).fetchone()
    if batch is None:
        abort(404)
    return batch


def _ok_report_or_404(db, bid, rid):
    report = db.execute(
        "SELECT * FROM reports WHERE id=? AND batch_id=?", (rid, bid)
    ).fetchone()
    if report is None or report["status"] != "ok":
        abort(404)
    return report


def _report_basis(db, report):
    """报告固定依据：所用 FASTA 版本解析出的 (chrom, ref_hash, ref_len, fasta_text)。"""
    fasta = db.execute(
        "SELECT content FROM inputs WHERE batch_id=? AND kind='fasta' AND version=?",
        (report["batch_id"], report["fasta_version"]),
    ).fetchone()
    if fasta is None:  # 理论上不会发生：报告与输入版本绑定
        abort(500)
    chrom, ref_hash, ref_len = annotate.ref_basis(fasta["content"])
    return chrom, ref_hash, ref_len, fasta["content"]


def _schema_columns(db, bid):
    row = db.execute(
        "SELECT columns FROM annotation_schemas WHERE batch_id=?", (bid,)
    ).fetchone()
    return json.loads(row["columns"]) if row else None


def _coverage(db, report, payload):
    """返回 (覆盖组下标集合, 位点 -> 注释 dict, 注释列)。缺口即未覆盖的标准位点。"""
    bid = report["batch_id"]
    chrom, ref_hash, _ref_len, _ = _report_basis(db, report)
    columns = _schema_columns(db, bid) or []
    rows = db.execute(
        "SELECT pos, ref, alt, data FROM annotations"
        " WHERE batch_id=? AND chrom=? AND ref_hash=?",
        (bid, chrom, ref_hash),
    ).fetchall()
    ann_by_key = {(r["pos"], r["ref"], r["alt"]): json.loads(r["data"]) for r in rows}
    covered = {i for i, g in enumerate(payload["groups"])
               if (g["pos"], g["ref"], g["alt"]) in ann_by_key}
    return covered, ann_by_key, columns


def _record_rejected(db, bid, rid, filename, rows_seen, error):
    """失败导入只写审计记录，不写任何注释数据（独立事务）。"""
    db.execute(
        "INSERT INTO annotation_imports(batch_id, report_id, filename, status,"
        " rows_seen, rows_inserted, rows_skipped, error, created_at)"
        " VALUES (?,?,?, 'rejected', ?, 0, 0, ?, datetime('now'))",
        (bid, rid, filename, rows_seen, error),
    )
    db.commit()


# ---------------------------------------------------------------- 页面

@bp.get("/batches/<int:bid>/reports/<int:rid>/annotations")
def annotation_page(bid, rid):
    db = get_db()
    batch = _batch_or_404(db, bid)
    report = _ok_report_or_404(db, bid, rid)
    payload = json.loads(report["payload"])
    covered, ann_by_key, columns = _coverage(db, report, payload)
    imports = db.execute(
        "SELECT * FROM annotation_imports WHERE batch_id=? ORDER BY id DESC LIMIT 30",
        (bid,),
    ).fetchall()
    fv, cv = (report["fasta_version"], report["csv_version"])
    cur_rows = db.execute(
        "SELECT kind, MAX(version) AS v FROM inputs WHERE batch_id=? GROUP BY kind",
        (bid,),
    ).fetchall()
    cur_v = {r["kind"]: r["v"] for r in cur_rows}
    is_current = (fv == cur_v.get("fasta") and cv == cur_v.get("csv"))
    return render_template(
        "annotations.html", batch=batch, report=report, payload=payload,
        covered=covered, ann_by_key=ann_by_key, columns=columns,
        imports=imports, is_current=is_current,
        fv=cur_v.get("fasta"), cv=cur_v.get("csv"),
    )


@bp.post("/batches/<int:bid>/reports/<int:rid>/annotations/import")
def import_annotation(bid, rid):
    db = get_db()
    _batch_or_404(db, bid)
    report = _ok_report_or_404(db, bid, rid)
    payload = json.loads(report["payload"])
    chrom, ref_hash, ref_len, _ = _report_basis(db, report)

    f = request.files.get("file")
    filename = f.filename if f is not None else ""
    if f is None or not filename:
        flash("请选择要导入的注释 CSV 文件")
        return redirect(url_for("annotations.annotation_page", bid=bid, rid=rid))
    try:
        text = f.read().decode("utf-8-sig")
    except UnicodeDecodeError:
        _record_rejected(db, bid, rid, filename, 0, "文件必须是 UTF-8 编码的文本")
        flash("导入被整次拒绝：文件必须是 UTF-8 编码的文本")
        return redirect(url_for("annotations.annotation_page", bid=bid, rid=rid))

    # 解析与业务校验在事务外完成：失败只留 rejected 审计，绝不触碰 annotations
    rows = []
    try:
        columns, rows = annotate.parse_annotation_csv(text)
        site_idx = annotate.standard_site_index(payload)
        annotate.check_against_sites(rows, site_idx, chrom)
        schema = _schema_columns(db, bid)
        if schema is not None and columns != schema:
            raise ValidationError(
                "注释列与本批次已固定的列不一致（整次拒绝）："
                f"本次为 {', '.join(columns)}；已固定为 {', '.join(schema)}"
            )
    except ValidationError as exc:
        _record_rejected(db, bid, rid, filename, len(rows), str(exc))
        flash(f"导入被整次拒绝：{exc}")
        return redirect(url_for("annotations.annotation_page", bid=bid, rid=rid))

    # 单事务写入：同值重传跳过、异值重传整次回滚、失败不留半次导入
    inserted = skipped = 0
    try:
        db.execute("BEGIN IMMEDIATE")
        if schema is None:
            db.execute(
                "INSERT INTO annotation_schemas(batch_id, columns, created_at)"
                " VALUES (?,?, datetime('now'))",
                (bid, json.dumps(columns, ensure_ascii=False)),
            )
        cur = db.execute(
            "INSERT INTO annotation_imports(batch_id, report_id, filename, status,"
            " rows_seen, rows_inserted, rows_skipped, error, created_at)"
            " VALUES (?,?,?, 'accepted', ?, 0, 0, NULL, datetime('now'))",
            (bid, rid, filename, len(rows)),
        )
        import_id = cur.lastrowid
        for r in rows:
            key = (bid, chrom, ref_hash, r["pos"], r["ref"], r["alt"])
            existing = db.execute(
                "SELECT data FROM annotations"
                " WHERE batch_id=? AND chrom=? AND ref_hash=? AND pos=? AND ref=? AND alt=?",
                key,
            ).fetchone()
            data_json = json.dumps(r["data"], ensure_ascii=False)
            if existing is not None:
                if existing["data"] != data_json:
                    raise ValidationError(
                        f"第 {r['line']} 行: 位点 {r['chrom']}:{r['pos']}:"
                        f"{r['ref']}>{r['alt']} 已存在不同注释值（同位点异值重传整次拒绝）"
                    )
                skipped += 1  # 同值重传：不重复入库
                continue
            db.execute(
                "INSERT INTO annotations(batch_id, chrom, ref_hash, ref_len,"
                " pos, ref, alt, data, import_id, created_at)"
                " VALUES (?,?,?,?,?,?,?,?,?, datetime('now'))",
                (bid, chrom, ref_hash, ref_len, r["pos"], r["ref"], r["alt"],
                 data_json, import_id),
            )
            inserted += 1
        db.execute(
            "UPDATE annotation_imports SET rows_inserted=?, rows_skipped=? WHERE id=?",
            (inserted, skipped, import_id),
        )
        db.commit()
    except Exception as exc:
        db.rollback()
        reason = str(exc) if isinstance(exc, ValidationError) else f"写入失败，已整体回滚: {exc}"
        _record_rejected(db, bid, rid, filename, len(rows), reason)
        flash(f"导入被整次拒绝：{reason}")
        return redirect(url_for("annotations.annotation_page", bid=bid, rid=rid))

    msg = f"导入完成：{len(rows)} 行，新入库 {inserted} 条"
    if skipped:
        msg += f"，同值重传跳过 {skipped} 条（未重复入库）"
    flash(msg)
    return redirect(url_for("annotations.annotation_page", bid=bid, rid=rid))


@bp.get("/batches/<int:bid>/reports/<int:rid>/annotations/download")
def download_annotation(bid, rid):
    """仅当报告全部标准位点都有注释时，才可下载按源 ID 展开的注释表。"""
    db = get_db()
    _batch_or_404(db, bid)
    report = _ok_report_or_404(db, bid, rid)
    payload = json.loads(report["payload"])
    covered, ann_by_key, columns = _coverage(db, report, payload)
    total = len(payload["groups"])
    if len(covered) != total:
        abort(409, description=f"仍有 {total - len(covered)} 个标准位点缺少注释，不能下载展开表")
    csv_text = annotate.expanded_annotation_csv(payload, ann_by_key, columns)
    return Response(
        csv_text,
        content_type="text/csv; charset=utf-8",
        headers={"Content-Disposition":
                 f"attachment; filename=batch{bid}_report{rid}_annotations.csv"},
    )
