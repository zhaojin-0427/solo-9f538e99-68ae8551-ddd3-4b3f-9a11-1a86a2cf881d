"""注释导入、回填复用、缺口与下载、并发交错与失败原子性测试。"""
import io
import json
import time

import pytest

from app import create_app
from app.annotate import (expanded_annotation_csv, parse_annotation_csv,
                          ref_basis, standard_site_index, check_against_sites)
from app.db import connect
from app.normalize import ValidationError
from app.worker import claim_job, process_job, run_pending_jobs

FASTA = ">chr1\n" + "A" * 20 + "\n"
CSV_OK = "ID,POS,REF,ALT\ni1,3,A,AA\ni2,9,A,AA\ns1,5,A,C\n"
# 标准化后两个标准位点：(1,A,AA) 组 ids=[i1,i2]，(5,A,C) 组 ids=[s1]
CSV_OK2 = "ID,POS,REF,ALT\nd1,10,AA,A\n"   # 更正后：标准位点 (1,AA,A)，与原位点不同
FASTA_V2 = ">chr1\n" + "A" * 30 + "\n"      # 染色体名相同、完整参考序列不同


@pytest.fixture()
def app(tmp_path):
    return create_app({
        "TESTING": True,
        "DATABASE": str(tmp_path / "t.db"),
        "WORKER_ENABLED": False,
        "COMPUTE_DELAY_SECONDS": 0,
    })


@pytest.fixture()
def client(app):
    return app.test_client()


def _upload(client, fasta=FASTA, csv_text=CSV_OK):
    return client.post("/batches", data={
        "name": "t",
        "fasta": (io.BytesIO(fasta.encode()), "ref.fa"),
        "csv": (io.BytesIO(csv_text.encode()), "vars.csv"),
    })


def _correct(client, bid, kind, text):
    return client.post(f"/batches/{bid}/inputs/{kind}", data={
        "file": (io.BytesIO(text.encode()), f"new.{kind}"),
    })


def _ann_upload(client, bid, rid, csv_text, filename="ann.csv"):
    return client.post(
        f"/batches/{bid}/reports/{rid}/annotations/import",
        data={"file": (io.BytesIO(csv_text.encode()), filename)},
        follow_redirects=True,
    )


def _report_payload(app, rid=1):
    conn = connect(app.config["DATABASE"])
    row = conn.execute("SELECT payload FROM reports WHERE id=?", (rid,)).fetchone()
    conn.close()
    return json.loads(row["payload"])


def _count_annotations(app):
    conn = connect(app.config["DATABASE"])
    n = conn.execute("SELECT COUNT(*) c FROM annotations").fetchone()["c"]
    conn.close()
    return n


def _imports(app):
    conn = connect(app.config["DATABASE"])
    rows = conn.execute(
        "SELECT * FROM annotation_imports ORDER BY id").fetchall()
    conn.close()
    return rows


FULL_ANN = (
    "CHROM,POS,REF,ALT,GENE,EFFECT\n"
    "chr1,1,A,AA,G1,ins_effect\n"
    "chr1,5,A,C,G2,snv_effect\n"
)


@pytest.fixture()
def ready(app, client):
    """建好批次并完成 v1 计算，返回 (bid, rid)。"""
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    return 1, 1


# ---------------------------------------------------------------- 纯函数

def test_parse_annotation_csv_basic():
    cols, rows = parse_annotation_csv(
        "chrom,pos,ref,alt,gene,score\nchr1,1,A,AA,G1,0.5\n")
    assert cols == ["gene", "score"]
    assert rows[0]["chrom"] == "chr1" and rows[0]["pos"] == 1
    assert rows[0]["data"] == {"gene": "G1", "score": "0.5"}


def test_parse_annotation_csv_requires_annotation_columns():
    with pytest.raises(ValidationError, match="至少需要一个注释列"):
        parse_annotation_csv("CHROM,POS,REF,ALT\nchr1,1,A,AA\n")


def test_parse_annotation_csv_rejects_duplicate_site():
    text = (
        "CHROM,POS,REF,ALT,G\n"
        "chr1,1,A,AA,G1\nchr1,1,A,AA,G2\n"
    )
    with pytest.raises(ValidationError, match="重复"):
        parse_annotation_csv(text)


def test_parse_annotation_csv_rejects_bad_pos_and_bases():
    with pytest.raises(ValidationError, match="正整数"):
        parse_annotation_csv("CHROM,POS,REF,ALT,G\nchr1,0,A,AA,G1\n")
    with pytest.raises(ValidationError, match="非法字符"):
        parse_annotation_csv("CHROM,POS,REF,ALT,G\nchr1,1,A,AX,G1\n")


def test_check_against_sites_unknown_and_wrong_chrom():
    payload = _payload()
    idx = standard_site_index(payload)
    good = [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "AA", "line": 2, "data": {}}]
    check_against_sites(good, idx, "chr1")  # 不抛异常
    bad_site = [{"chrom": "chr1", "pos": 2, "ref": "A", "alt": "AA", "line": 2, "data": {}}]
    with pytest.raises(ValidationError, match="未知位点"):
        check_against_sites(bad_site, idx, "chr1")
    bad_chrom = [{"chrom": "chr2", "pos": 1, "ref": "A", "alt": "AA", "line": 2, "data": {}}]
    with pytest.raises(ValidationError, match="染色体名"):
        check_against_sites(bad_chrom, idx, "chr1")


def test_ref_basis_changes_with_sequence_not_name():
    c1, h1, n1 = ref_basis(FASTA)
    c2, h2, n2 = ref_basis(FASTA_V2)
    assert c1 == c2 == "chr1"          # 染色体名相同
    assert h1 != h2 and n1 != n2       # 完整参考序列不同 → 依据不同


def test_expanded_csv_one_row_per_source_id():
    payload = _payload()
    ann = {(1, "A", "AA"): {"GENE": "G1", "EFFECT": "e1"},
           (5, "A", "C"): {"GENE": "G2", "EFFECT": "e2"}}
    text = expanded_annotation_csv(payload, ann, ["GENE", "EFFECT"])
    lines = text.strip().splitlines()
    assert lines[0] == "ID,CHROM,POS,REF,ALT,GENE,EFFECT"
    # 三个源 ID 各一行；i1/i2 同组共享同一条注释
    data = {ln.split(",")[0]: ln for ln in lines[1:]}
    assert set(data) == {"i1", "i2", "s1"}
    # 同组的 i1/i2 除 ID 外内容完全一致，共享同一条注释
    assert data["i1"][len("i1"):] == data["i2"][len("i2"):] == ",chr1,1,A,AA,G1,e1"
    assert data["s1"] == "s1,chr1,5,A,C,G2,e2"


def _payload():
    """复用 worker 计算出的标准 payload（等价归并后的两组）。"""
    from app.normalize import compute_report
    return compute_report(FASTA, CSV_OK)


# ---------------------------------------------------------------- 导入基本流程

def test_import_full_coverage_enables_download(app, client, ready):
    bid, rid = ready
    page = _ann_upload(client, bid, rid, FULL_ANN).data.decode()
    assert "全部位点已覆盖" in page and "已注释" in page
    assert _count_annotations(app) == 2

    resp = client.get(f"/batches/{bid}/reports/{rid}/annotations/download")
    assert resp.status_code == 200
    body = resp.data.decode()
    lines = body.strip().splitlines()
    assert lines[0] == "ID,CHROM,POS,REF,ALT,GENE,EFFECT"
    assert len(lines) == 4                     # 表头 + 3 个源 ID
    assert "i1,chr1,1,A,AA,G1,ins_effect" in body
    assert "i2,chr1,1,A,AA,G1,ins_effect" in body

    imps = _imports(app)
    assert len(imps) == 1 and imps[0]["status"] == "accepted"
    assert imps[0]["rows_seen"] == 2 and imps[0]["rows_inserted"] == 2


def test_partial_coverage_blocks_download_and_shows_gaps(app, client, ready):
    bid, rid = ready
    partial = "CHROM,POS,REF,ALT,GENE\nchr1,5,A,C,G2\n"
    page = _ann_upload(client, bid, rid, partial).data.decode()
    assert "缺口" in page and "1 / 2" in page

    resp = client.get(f"/batches/{bid}/reports/{rid}/annotations/download")
    assert resp.status_code == 409
    # 补齐后可下载，且不会漏掉先到的注释
    rest = "CHROM,POS,REF,ALT,GENE\nchr1,1,A,AA,G1\n"
    _ann_upload(client, bid, rid, rest)
    resp = client.get(f"/batches/{bid}/reports/{rid}/annotations/download")
    assert resp.status_code == 200
    body = resp.data.decode()
    assert "i1,chr1,1,A,AA,G1" in body and "s1,chr1,5,A,C,G2" in body


def test_unknown_site_rejects_whole_import_without_partial(app, client, ready):
    bid, rid = ready
    text = (
        "CHROM,POS,REF,ALT,GENE\n"
        "chr1,5,A,C,G2\n"          # 已知
        "chr1,9,A,AA,GX\n"        # 未知（标准位点是 1:A:AA，不是 9）
    )
    page = _ann_upload(client, bid, rid, text).data.decode()
    assert "整次拒绝" in page and "未知位点" in page
    assert _count_annotations(app) == 0        # 不留半次导入
    imps = _imports(app)
    assert len(imps) == 1 and imps[0]["status"] == "rejected"
    assert imps[0]["rows_inserted"] == 0 and imps[0]["error"]


def test_duplicate_site_in_file_rejects_whole_import(app, client, ready):
    bid, rid = ready
    text = (
        "CHROM,POS,REF,ALT,GENE\n"
        "chr1,1,A,AA,G1\nchr1,1,A,AA,G1b\n"
    )
    page = _ann_upload(client, bid, rid, text).data.decode()
    assert "整次拒绝" in page and "重复" in page
    assert _count_annotations(app) == 0


def test_wrong_chrom_rejects(app, client, ready):
    bid, rid = ready
    text = "CHROM,POS,REF,ALT,GENE\nchr2,1,A,AA,G1\n"
    page = _ann_upload(client, bid, rid, text).data.decode()
    assert "整次拒绝" in page and "染色体名" in page
    assert _count_annotations(app) == 0


def test_same_value_reupload_is_idempotent(app, client, ready):
    bid, rid = ready
    _ann_upload(client, bid, rid, FULL_ANN)
    page = _ann_upload(client, bid, rid, FULL_ANN).data.decode()
    assert "同值重传跳过 2 条" in page
    assert _count_annotations(app) == 2        # 不重复入库
    imps = _imports(app)
    assert imps[-1]["status"] == "accepted"
    assert imps[-1]["rows_inserted"] == 0 and imps[-1]["rows_skipped"] == 2


def test_different_value_same_site_rejects_everything(app, client, ready):
    bid, rid = ready
    _ann_upload(client, bid, rid, FULL_ANN)
    changed = (
        "CHROM,POS,REF,ALT,GENE,EFFECT\n"
        "chr1,1,A,AA,G1,ins_effect\n"
        "chr1,5,A,C,G2,DIFFERENT\n"   # 同位点异值
    )
    page = _ann_upload(client, bid, rid, changed).data.decode()
    assert "整次拒绝" in page and "不同注释值" in page
    assert _count_annotations(app) == 2        # 原有数据不动，本次零写入
    # 原值仍可下载
    body = client.get(f"/batches/{bid}/reports/{rid}/annotations/download").data.decode()
    assert "snv_effect" in body and "DIFFERENT" not in body


def test_schema_columns_fixed_per_batch(app, client, ready):
    bid, rid = ready
    _ann_upload(client, bid, rid, "CHROM,POS,REF,ALT,GENE\nchr1,5,A,C,G2\n")
    other_cols = "CHROM,POS,REF,ALT,SCORE\nchr1,1,A,AA,9\n"
    page = _ann_upload(client, bid, rid, other_cols).data.decode()
    assert "整次拒绝" in page and "注释列" in page
    assert _count_annotations(app) == 1


def test_failed_report_has_no_annotation_page(app, client):
    _upload(client, csv_text="ID,POS,REF,ALT\nbad,5,C,G\n")
    run_pending_jobs(app.config["DATABASE"])
    resp = client.get("/batches/1/reports/1/annotations")
    assert resp.status_code == 404


# ---------------------------------------------------------------- 版本复用与缺口

def test_new_report_reuses_only_when_chrom_ref_and_site_match(app, client, ready):
    bid, rid = ready
    _ann_upload(client, bid, rid, FULL_ANN)

    # CSV 更正：染色体名/参考序列相同，但标准位点变为 (1,AA,A) → 全部缺口
    _correct(client, bid, "csv", CSV_OK2)
    run_pending_jobs(app.config["DATABASE"])
    page = client.get("/batches/1/reports/2/annotations").data.decode()
    assert "0 / 1" in page and "缺口" in page
    assert client.get("/batches/1/reports/2/annotations/download").status_code == 409
    # 注释数据未被删除，旧报告仍完整可查可下载（固定原依据，标过期）
    old = client.get("/batches/1/reports/1/annotations").data.decode()
    assert "已过期" in old and "全部位点已覆盖" in old
    assert client.get("/batches/1/reports/1/annotations/download").status_code == 200


def test_fasta_correction_same_chrom_but_different_ref_creates_gaps(app, client, ready):
    bid, rid = ready
    _ann_upload(client, bid, rid, FULL_ANN)
    # 参考序列更正：染色体名仍为 chr1，但完整序列哈希变化 → 旧注释不可复用
    _correct(client, bid, "fasta", FASTA_V2)
    run_pending_jobs(app.config["DATABASE"])
    page = client.get("/batches/1/reports/2/annotations").data.decode()
    assert "0 / 2" in page and "缺口" in page
    # 旧报告固定原依据，仍可下载
    assert client.get("/batches/1/reports/1/annotations/download").status_code == 200


def test_annotation_imported_against_old_report_reused_by_matching_new_report(
        app, client, ready):
    """先更正产生新报告，再给旧报告导入；只要依据与位点相同，新报告即时复用。"""
    bid, rid = ready
    # 一次 CSV 更正后再更正回与 v1 等价的内容：报告 v3 位点与 v1 相同
    _correct(client, bid, "csv", CSV_OK2)
    run_pending_jobs(app.config["DATABASE"])
    _correct(client, bid, "csv", CSV_OK)
    run_pending_jobs(app.config["DATABASE"])
    # 向已过期的报告 #1 导入注释
    page = _ann_upload(client, bid, 1, FULL_ANN).data.decode()
    assert "已过期" in page
    # 当前报告 #3 同染色体、同完整参考、同标准位点 → 自动复用，立即可下载
    cur = client.get("/batches/1/reports/3/annotations").data.decode()
    assert "全部位点已覆盖" in cur
    resp = client.get("/batches/1/reports/3/annotations/download")
    assert resp.status_code == 200 and "i1,chr1,1,A,AA,G1,ins_effect" in resp.data.decode()


def test_chrom_rename_breaks_reuse(app, client, ready):
    bid, rid = ready
    _ann_upload(client, bid, rid, FULL_ANN)
    _correct(client, bid, "fasta", ">chr1_renamed\n" + "A" * 20 + "\n")
    run_pending_jobs(app.config["DATABASE"])
    page = client.get("/batches/1/reports/2/annotations").data.decode()
    assert "0 / 2" in page  # 序列相同但染色体名不同 → 不复用


# ---------------------------------------------------------------- 并发交错

def test_import_during_background_recompute_never_treats_old_basis_as_current(
        tmp_path):
    """真实 worker 线程：重算进行中导入注释，新报告缺口/复用按依据即时判定。"""
    db_path = str(tmp_path / "live.db")
    app = create_app({
        "TESTING": True, "DATABASE": db_path, "WORKER_ENABLED": True,
        "COMPUTE_DELAY_SECONDS": 1.0, "WORKER_POLL_SECONDS": 0.05,
    })
    client = app.test_client()
    _upload(client)
    time.sleep(0.3)                                    # v1 计算中
    # 计算期间提交 CSV 更正（位点将变化）
    _correct(client, 1, "csv", CSV_OK2)
    deadline = time.time() + 10
    while time.time() < deadline:
        conn = connect(db_path)
        n = conn.execute("SELECT COUNT(*) c FROM reports").fetchone()["c"]
        conn.close()
        if n >= 2:
            break
        time.sleep(0.1)
    # 给旧报告导入（旧依据），新报告位点不同 → 必须显示缺口，不能当成当前完整结果
    page = _ann_upload(client, 1, 1,
                       "CHROM,POS,REF,ALT,GENE\nchr1,1,A,AA,G1\nchr1,5,A,C,G2\n"
                       ).data.decode()
    assert "已过期" in page
    new_page = client.get("/batches/1/reports/2/annotations").data.decode()
    assert "0 / 1" in new_page and "缺口" in new_page
    assert client.get("/batches/1/reports/2/annotations/download").status_code == 409
    # 旧报告依据固定，仍完整可下载
    assert client.get("/batches/1/reports/1/annotations/download").status_code == 200


def test_interleaved_imports_and_corrections_lose_no_reusable_annotation(app, client, ready):
    """导入与更正任意顺序交错：满足依据的历史注释一条都不漏。"""
    bid, rid = ready
    # 先导半个（仅 SNV 位点）
    _ann_upload(client, bid, rid, "CHROM,POS,REF,ALT,GENE\nchr1,5,A,C,G2\n")
    # 更正又改回（报告 v3 与 v1 同依据同位点）
    _correct(client, bid, "csv", CSV_OK2)
    run_pending_jobs(app.config["DATABASE"])
    _correct(client, bid, "csv", CSV_OK)
    run_pending_jobs(app.config["DATABASE"])
    # 再补另一半
    page = _ann_upload(client, bid, 3,
                       "CHROM,POS,REF,ALT,GENE\nchr1,1,A,AA,G1\n").data.decode()
    assert "全部位点已覆盖" in page
    body = client.get("/batches/1/reports/3/annotations/download").data.decode()
    assert "G1" in body and "G2" in body              # 两批注释都被复用


def test_failed_import_writes_no_half_data_and_keeps_audit(app, client, ready):
    bid, rid = ready
    _ann_upload(client, bid, rid, "CHROM,POS,REF,ALT,GENE\nchr1,5,A,C,G2\n")
    # 混合：一行合法、一行异值 → 整次回滚，合法行也不得入库
    text = (
        "CHROM,POS,REF,ALT,GENE\n"
        "chr1,1,A,AA,G1\n"
        "chr1,5,A,C,CHANGED\n"
    )
    _ann_upload(client, bid, rid, text)
    conn = connect(app.config["DATABASE"])
    n = conn.execute("SELECT COUNT(*) c FROM annotations").fetchone()["c"]
    imps = conn.execute("SELECT status, rows_inserted FROM annotation_imports").fetchall()
    conn.close()
    assert n == 1                                     # 只有第一次的那条
    assert [r["status"] for r in imps] == ["accepted", "rejected"]
    assert imps[1]["rows_inserted"] == 0


def test_batch_page_shows_annotation_status(app, client, ready):
    bid, rid = ready
    page = client.get(f"/batches/{bid}").data.decode()
    assert "注释覆盖" in page and "0 / 2" in page
    _ann_upload(client, bid, rid, FULL_ANN)
    page = client.get(f"/batches/{bid}").data.decode()
    assert "已齐全" in page and "下载按源 ID 展开的注释表" in page


def test_non_utf8_file_rejected_with_audit(app, client, ready):
    bid, rid = ready
    resp = client.post(
        f"/batches/{bid}/reports/{rid}/annotations/import",
        data={"file": (io.BytesIO(b"CHROM,POS,REF,ALT,G\nchr1,\xff\xfe"), "bad.csv")},
        follow_redirects=True,
    )
    assert "整次拒绝" in resp.data.decode() and "UTF-8" in resp.data.decode()
    assert _count_annotations(app) == 0
    imps = _imports(app)
    assert len(imps) == 1 and imps[0]["status"] == "rejected"


def test_missing_file_does_not_create_audit(app, client, ready):
    bid, rid = ready
    resp = client.post(
        f"/batches/{bid}/reports/{rid}/annotations/import",
        data={}, follow_redirects=True,
    )
    assert resp.status_code == 200 and "请选择" in resp.data.decode()
    assert len(_imports(app)) == 0 and _count_annotations(app) == 0


def test_annotation_page_unknown_batch_or_report_404(app, client, ready):
    assert client.get("/batches/999/reports/1/annotations").status_code == 404
    assert client.get("/batches/1/reports/999/annotations").status_code == 404


def test_concurrent_imports_no_duplicate_rows(app, client, ready):
    """两个并发请求写同一位点：UNIQUE 约束 + 单事务保证不产生重复注释。"""
    import threading

    bid, rid = ready
    text = "CHROM,POS,REF,ALT,GENE\nchr1,5,A,C,G2\n"
    results = []

    def fire():
        results.append(_ann_upload(client, bid, rid, text).status_code)

    threads = [threading.Thread(target=fire) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert results == [200] * 4
    conn = connect(app.config["DATABASE"])
    n = conn.execute(
        "SELECT COUNT(*) c FROM annotations WHERE pos=5 AND ref='A' AND alt='C'"
    ).fetchone()["c"]
    accepted = conn.execute(
        "SELECT COUNT(*) c FROM annotation_imports WHERE status='accepted'"
    ).fetchone()["c"]
    rejected = conn.execute(
        "SELECT COUNT(*) c FROM annotation_imports WHERE status='rejected'"
    ).fetchone()["c"]
    conn.close()
    assert n == 1                    # 无论如何只有一条注释
    assert accepted + rejected == 4  # 每次导入都有审计，无半次、无丢失
