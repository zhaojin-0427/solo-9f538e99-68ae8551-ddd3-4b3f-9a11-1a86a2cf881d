"""注释 CSV 解析、标准位点匹配与按源 ID 展开（纯函数模块）。

不依赖 Flask / SQLite，便于单独测试。非法输入抛 normalize.ValidationError，
消息面向用户展示；调用方保证整次导入失败、不留半次数据。

匹配键：(CHROM, POS, REF, ALT)，其中 POS/REF/ALT 取自标准化报告的标准位点。
"""
from __future__ import annotations

import csv
import hashlib
import io

from .normalize import BASES, ValidationError, parse_fasta

MAX_ANN_ROWS = 1000
KEY_COLS = ("chrom", "pos", "ref", "alt")


# ---------------------------------------------------------------- 解析

def parse_annotation_csv(text: str):
    """解析注释 CSV。

    返回 (ann_columns, rows)：
    - ann_columns: 四个键列之外的注释列名（保持文件中的顺序，批次内固定）；
    - rows: [{"chrom", "pos", "ref", "alt", "data": {列: 值}, "line": 文件行号}]。
    """
    try:
        raw = [r for r in csv.reader(io.StringIO(text))]
    except csv.Error as exc:
        raise ValidationError(f"CSV 解析失败: {exc}") from exc
    raw = [r for r in raw if any(cell.strip() for cell in r)]
    if not raw:
        raise ValidationError("注释 CSV 文件为空")

    header = [c.strip() for c in raw[0]]
    if any(not h for h in header):
        raise ValidationError("注释 CSV 表头存在空列名")
    if len(set(header)) != len(header):
        dup = sorted({h for h in header if header.count(h) > 1})
        raise ValidationError("注释 CSV 表头列名重复: " + ", ".join(dup))
    missing = [c.upper() for c in KEY_COLS if c not in [h.lower() for h in header]]
    if missing:
        raise ValidationError("注释 CSV 缺少必需列: " + ", ".join(missing))

    lower = [h.lower() for h in header]
    idx = {c: lower.index(c) for c in KEY_COLS}
    ann_columns = [h for h, low in zip(header, lower) if low not in KEY_COLS]
    if not ann_columns:
        raise ValidationError(
            "注释 CSV 只有 CHROM/POS/REF/ALT 四个键列，至少需要一个注释列"
        )

    data_rows = raw[1:]
    if not data_rows:
        raise ValidationError("注释 CSV 没有数据行")
    if len(data_rows) > MAX_ANN_ROWS:
        raise ValidationError(f"注释 CSV 数据行数 {len(data_rows)} 超过上限 {MAX_ANN_ROWS}")

    out, seen = [], set()
    for n, r in enumerate(data_rows, start=2):  # n = 含表头的文件行号
        def cell(name, r=r):
            i = idx[name]
            return r[i].strip() if i < len(r) else ""

        chrom = cell("chrom")
        if not chrom:
            raise ValidationError(f"第 {n} 行: CHROM 为空")
        pos_s = cell("pos")
        if not pos_s.isdigit() or int(pos_s) < 1:
            raise ValidationError(f"第 {n} 行: POS '{pos_s}' 不是从 1 开始的正整数")
        pos = int(pos_s)
        ref, alt = cell("ref").upper(), cell("alt").upper()
        if not ref or not alt:
            raise ValidationError(f"第 {n} 行: REF/ALT 不能为空")
        for label, v in (("REF", ref), ("ALT", alt)):
            bad = sorted(set(v) - BASES)
            if bad:
                raise ValidationError(f"第 {n} 行: {label} 含非法字符 {''.join(bad)}")

        key = (chrom, pos, ref, alt)
        if key in seen:
            raise ValidationError(
                f"第 {n} 行: 位点 {chrom}:{pos}:{ref}>{alt} 在文件内重复（重复位点整次拒绝）"
            )
        seen.add(key)

        data = {}
        for col in ann_columns:
            i = header.index(col)
            data[col] = r[i].strip() if i < len(r) else ""
        out.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "data": data, "line": n,
        })
    return ann_columns, out


# ---------------------------------------------------------------- 位点匹配

def standard_site_index(payload: dict) -> dict:
    """从标准化报告 payload 构造标准位点索引 {(pos, ref, alt): group}。"""
    return {(g["pos"], g["ref"], g["alt"]): g for g in payload["groups"]}


def check_against_sites(rows: list, site_idx: dict, chrom: str) -> None:
    """逐行核对 CHROM 与标准位点；未知位点整次拒绝。"""
    for r in rows:
        if r["chrom"] != chrom:
            raise ValidationError(
                f"第 {r['line']} 行: CHROM '{r['chrom']}' 与报告染色体名 '{chrom}' 不一致"
            )
        key = (r["pos"], r["ref"], r["alt"])
        if key not in site_idx:
            raise ValidationError(
                f"第 {r['line']} 行: 位点 {r['chrom']}:{r['pos']}:"
                f"{r['ref']}>{r['alt']} 不在该报告的标准位点中（未知位点整次拒绝）"
            )


# ---------------------------------------------------------------- 参考依据

def ref_basis(fasta_text: str):
    """注释复用依据：(染色体名, 完整参考序列 SHA256, 参考长度)。"""
    chrom, seq = parse_fasta(fasta_text)
    digest = hashlib.sha256(seq.encode("ascii")).hexdigest()
    return chrom, digest, len(seq)


# ---------------------------------------------------------------- 展开下载

def expanded_annotation_csv(payload: dict, ann_by_key: dict, columns: list) -> str:
    """按源 ID 展开注释表：每个源 ID 一行（等价归并组共享同一条注释）。"""
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["ID", "CHROM", "POS", "REF", "ALT", *columns])
    for src in payload["rows"]:
        group = payload["groups"][src["group"]]
        key = (group["pos"], group["ref"], group["alt"])
        ann = ann_by_key.get(key, {})
        writer.writerow([
            src["id"], payload["chrom"],
            group["pos"], group["ref"], group["alt"],
            *(ann.get(c, "") for c in columns),
        ])
    return buf.getvalue()
