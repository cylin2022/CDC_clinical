#!/usr/bin/env python3
"""Offline Clinical Nanopore Pathogen Platform — Web GUI & FIFO Job Queue Server."""

import base64
import csv
import datetime as dt
import glob
import html
import io
import json
import math
import mimetypes
import os
import posixpath
import queue
import re
import shutil
import sqlite3
import subprocess
import sys
import threading
import time
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT_DIR, "scripts"))
from resource_optimizer import compute_execution_plan, get_system_resources, warm_kraken2_cache  # noqa: E402

INPUT_DIR = os.environ.get("CDC_INPUT_DIR", os.path.join(ROOT_DIR, "input_samples"))
OUTPUT_DIR = os.environ.get("CDC_OUTPUT_DIR", os.path.join(ROOT_DIR, "output_results"))
DB_DIR = os.environ.get("CDC_DB_DIR", os.path.join(ROOT_DIR, "databases"))
REF_DIR = os.environ.get("CDC_REF_DIR", os.path.join(ROOT_DIR, "references"))
PORT = int(os.environ.get("CDC_PORT", "8080"))

JOBS_DB = os.path.join(OUTPUT_DIR, "jobs.db")
DB_LOCK = threading.Lock()
ACTIVE_PROCS = {}
PROC_LOCK = threading.Lock()


def init_db():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(INPUT_DIR, exist_ok=True)
    with DB_LOCK:
        conn = sqlite3.connect(JOBS_DB)
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sample_id TEXT NOT NULL,
                dna_path TEXT,
                cdna_path TEXT,
                specimen_type TEXT DEFAULT 'unknown',
                threads INTEGER DEFAULT 16,
                priority INTEGER DEFAULT 100,
                status TEXT NOT NULL DEFAULT 'queued',
                stage TEXT DEFAULT 'Queued',
                progress_pct INTEGER DEFAULT 0,
                outdir TEXT NOT NULL,
                primary_finding TEXT DEFAULT '',
                error_msg TEXT DEFAULT '',
                created_at TEXT NOT NULL,
                started_at TEXT DEFAULT '',
                finished_at TEXT DEFAULT ''
            )
            """
        )
        # Seed pre-existing sample_001 historical analysis if present
        cur = conn.execute("SELECT id FROM jobs WHERE sample_id = 'sample_001'").fetchone()
        val_dir = os.path.join(ROOT_DIR, "analysis", "pathogen_validation")
        if not cur and os.path.exists(os.path.join(val_dir, "read_disposition_sankey.html")):
            conn.execute(
                """
                INSERT INTO jobs (
                    sample_id, dna_path, cdna_path, specimen_type, threads, priority,
                    status, stage, progress_pct, outdir, primary_finding,
                    created_at, started_at, finished_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    "sample_001",
                    os.path.join(ROOT_DIR, "fastq_pass"),
                    os.path.join(ROOT_DIR, "fastq_pass"),
                    "clinical (DNA + cDNA)",
                    38,
                    1,
                    "completed",
                    "Completed",
                    100,
                    val_dir,
                    "Influenza A (H3N2) [98.45% cov, 78.1x] + Candida albicans [DNA 4.41%]",
                    "2026-11-12T21:00:00+0800",
                    "2026-11-12T21:00:00+0800",
                    "2026-11-13T12:12:04+0800",
                ),
            )
        conn.commit()
        conn.close()


def dict_factory(cursor, row):
    return {col[0]: row[idx] for idx, col in enumerate(cursor.description)}


def get_jobs():
    with DB_LOCK:
        conn = sqlite3.connect(JOBS_DB)
        conn.row_factory = dict_factory
        rows = conn.execute(
            "SELECT * FROM jobs ORDER BY CASE status WHEN 'running' THEN 0 WHEN 'queued' THEN 1 ELSE 2 END, priority ASC, id ASC"
        ).fetchall()
        conn.close()

    for r in rows:
        prog_file = os.path.join(r["outdir"], "progress.json")
        if r["status"] == "running" and os.path.exists(prog_file):
            try:
                with open(prog_file) as f:
                    p = json.load(f)
                r["stage"] = p.get("stage", r["stage"])
                r["progress_pct"] = p.get("pct", r["progress_pct"])
            except Exception:
                pass
    return rows


def dir_size_human(path):
    if not os.path.exists(path):
        return None
    if os.path.isfile(path):
        b = os.path.getsize(path)
    else:
        b = 0
        for dp, _, fns in os.walk(path):
            for fn in fns:
                fp = os.path.join(dp, fn)
                if not os.path.islink(fp):
                    b += os.path.getsize(fp)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if b < 1024.0:
            return f"{b:.1f} {unit}"
        b /= 1024.0
    return f"{b:.1f} PB"


def check_offline_env():
    k2_path = os.path.join(DB_DIR, "kraken2", "k2_pluspf_20261022", "hash.k2d")
    hm_dna = os.path.join(REF_DIR, "human-masked.mmi")
    hm_cdna = os.path.join(REF_DIR, "human-masked.splice.mmi")
    dmnd_db = os.path.join(DB_DIR, "viral_refseq_prot", "viral_refseq_prot.dmnd")
    panel_dir = os.path.join(REF_DIR, "pathogen_panel")
    panel_cnt = len(glob.glob(os.path.join(panel_dir, "*.fna.gz")))
    return {
        "offline_mode": True,
        "input_dir": INPUT_DIR,
        "output_dir": OUTPUT_DIR,
        "databases": [
            {
                "name": "Kraken2 PlusPF (2026-10-22)",
                "ready": os.path.exists(k2_path),
                "size": dir_size_human(os.path.dirname(k2_path)) or "Missing",
                "description": "k-mer 宏基因組快速物種分類主庫（涵蓋古菌、細菌、病毒、質體、人類、UniVec_Core、原蟲與真菌），支援 memory-mapping 記憶體映射加速。",
            },
            {
                "name": "Hostile Human Masked Indices (DNA + cDNA)",
                "ready": os.path.exists(hm_dna) and os.path.exists(hm_cdna),
                "size": dir_size_human(REF_DIR) or "Missing",
                "description": "T2T-CHM13v2.0 + rCRS 粒線體 + IPD-IMGT/HLA 人類去宿主索引（已遮蔽細菌/病毒/噬菌體保守區避免誤刪微生物）；含 DNA map-ont 與 cDNA splice 雙索引。",
            },
            {
                "name": "NCBI Viral RefSeq Protein (DIAMOND)",
                "ready": os.path.exists(dmnd_db),
                "size": dir_size_human(os.path.dirname(dmnd_db)) or "Missing",
                "description": "NCBI 病毒參考蛋白庫（DIAMOND blastx），用於對未分類序列進行蛋白層級遠緣同源搜尋，篩檢潛在未知或高度變異病毒。",
            },
            {
                "name": f"Offline Pathogen Panel ({panel_cnt} genomes)",
                "ready": panel_cnt > 0,
                "size": dir_size_human(panel_dir) or "Missing",
                "description": "離線病原體競爭性比對面板（收錄病毒、細菌、真菌及常見試劑污染代表株之完整參考基因組），用於回貼確認覆蓋廣度、深度與 MAPQ 信心值。",
                "has_catalog": os.path.exists(os.path.join(panel_dir, "PANEL_CATALOG.tsv")),
            },
        ],
    }


def scan_input_samples(scan_dir=None):
    """Scan scan_dir (or INPUT_DIR) for <sample_id>/{dna,cdna}/*.fastq.gz folders."""
    target_dir = scan_dir if scan_dir and os.path.isdir(scan_dir) else INPUT_DIR
    samples = []
    if not os.path.isdir(target_dir):
        return target_dir, samples
    for name in sorted(os.listdir(target_dir)):
        sdir = os.path.join(target_dir, name)
        if not os.path.isdir(sdir) or name.startswith("."):
            continue
        dna_dir = ""
        cdna_dir = ""
        dna_files = []
        cdna_files = []
        for cand in ("dna", "DNA"):
            p = os.path.join(sdir, cand)
            if os.path.isdir(p):
                fqs = sorted(glob.glob(os.path.join(p, "*.f*q*")))
                if fqs:
                    dna_dir = p
                    dna_files = [os.path.basename(x) for x in fqs]
                    break
        for cand in ("cdna", "cDNA", "rna", "RNA"):
            p = os.path.join(sdir, cand)
            if os.path.isdir(p):
                fqs = sorted(glob.glob(os.path.join(p, "*.f*q*")))
                if fqs:
                    cdna_dir = p
                    cdna_files = [os.path.basename(x) for x in fqs]
                    break
        if not dna_files and not cdna_files:
            # Also list direct fastq files if folder itself is a barcode/run folder
            direct_fqs = sorted(glob.glob(os.path.join(sdir, "*.f*q*")))
            if direct_fqs:
                dna_dir = sdir
                dna_files = [os.path.basename(x) for x in direct_fqs]
        samples.append({
            "sample_id": name,
            "dna_path": dna_dir,
            "cdna_path": cdna_dir,
            "dna_files": dna_files,
            "cdna_files": cdna_files,
        })
    return target_dir, samples


def browse_directory(req_path):
    """List subdirectories and FASTQ files in req_path for the Web GUI directory picker."""
    target = os.path.abspath(req_path) if req_path else ROOT_DIR
    if not os.path.isdir(target):
        target = ROOT_DIR
    parent = os.path.dirname(target) if target != "/" else "/"
    dirs = []
    files = []
    try:
        for entry in sorted(os.scandir(target), key=lambda e: e.name.lower()):
            if entry.name.startswith("."):
                continue
            if entry.is_dir(follow_symlinks=True):
                dirs.append({"name": entry.name, "path": entry.path})
            elif entry.is_file(follow_symlinks=True) and (
                entry.name.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))
            ):
                sz = entry.stat().st_size
                for u in ("B", "KB", "MB", "GB", "TB"):
                    if sz < 1024.0:
                        sz_str = f"{sz:.1f} {u}"
                        break
                    sz /= 1024.0
                else:
                    sz_str = f"{sz:.1f} PB"
                files.append({"name": entry.name, "path": entry.path, "size": sz_str})
    except PermissionError:
        pass
    quick_roots = [
        {"label": "工作目錄", "path": ROOT_DIR},
        {"label": "input_samples", "path": INPUT_DIR},
    ]
    fq_pass = os.path.join(ROOT_DIR, "fastq_pass")
    if os.path.isdir(fq_pass):
        quick_roots.append({"label": "fastq_pass", "path": fq_pass})
    if os.path.isdir("/data"):
        quick_roots.append({"label": "/data", "path": "/data"})
    quick_roots.append({"label": "/", "path": "/"})
    return {
        "current_path": target,
        "parent_path": parent,
        "quick_roots": quick_roots,
        "dirs": dirs,
        "files": files,
    }


def _run_job_worker(job, threads):
    runner = os.path.join(ROOT_DIR, "scripts", "run_offline_sample.py")
    cmd = [
        sys.executable, runner,
        "--sample-id", job["sample_id"],
        "--dna-path", job["dna_path"] or "",
        "--cdna-path", job["cdna_path"] or "",
        "--specimen-type", job["specimen_type"] or "unknown",
        "--threads", str(threads),
        "--db-dir", DB_DIR,
        "--ref-dir", REF_DIR,
        "--outdir", job["outdir"],
    ]
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        with PROC_LOCK:
            ACTIVE_PROCS[job["id"]] = proc
        out_log, _ = proc.communicate()
        rc = proc.returncode
    except Exception as exc:
        out_log = str(exc)
        rc = -1
    finally:
        with PROC_LOCK:
            ACTIVE_PROCS.pop(job["id"], None)

    fin = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    primary_finding = ""
    sum_json = os.path.join(job["outdir"], "sample_summary.json")
    if rc == 0 and os.path.exists(sum_json):
        try:
            with open(sum_json) as sf:
                primary_finding = json.load(sf).get("primary_finding", "")
        except Exception:
            pass
    with DB_LOCK:
        conn = sqlite3.connect(JOBS_DB)
        if rc == 0:
            conn.execute(
                "UPDATE jobs SET status='completed', stage='Completed', progress_pct=100, primary_finding=?, finished_at=? WHERE id=?",
                (primary_finding, fin, job["id"]),
            )
        else:
            conn.execute(
                "UPDATE jobs SET status='failed', stage='Failed', error_msg=?, finished_at=? WHERE id=?",
                (out_log[-1000:] if out_log else f"Exit code {rc}", fin, job["id"]),
            )
        conn.commit()
        conn.close()


def worker_loop():
    while True:
        with PROC_LOCK:
            active_count = len(ACTIVE_PROCS)
        with DB_LOCK:
            conn = sqlite3.connect(JOBS_DB)
            conn.row_factory = dict_factory
            queued_count = conn.execute("SELECT COUNT(*) AS c FROM jobs WHERE status = 'queued'").fetchone()["c"]
            plan = compute_execution_plan(queued_count=queued_count, active_count=active_count)
            job = None
            if queued_count > 0 and active_count < plan["max_workers"]:
                job = conn.execute(
                    "SELECT * FROM jobs WHERE status = 'queued' ORDER BY priority ASC, id ASC LIMIT 1"
                ).fetchone()
                if job:
                    now = time.strftime("%Y-%m-%dT%H:%M:%S%z")
                    threads = plan["threads_per_job"]
                    conn.execute(
                        "UPDATE jobs SET status='running', stage='Starting', progress_pct=5, threads=?, started_at=? WHERE id=?",
                        (threads, now, job["id"]),
                    )
                    conn.commit()
            conn.close()

        if not job:
            time.sleep(2)
            continue

        threading.Thread(target=_run_job_worker, args=(job, threads), daemon=True).start()
        time.sleep(0.3)


def read_tsv(path, limit=50):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))[:limit]


def _is_positive_finding(pf):
    pf = (pf or "").strip()
    return bool(pf) and not pf.lower().startswith("no major pathogen")


def _extract_pathogen_names(pf):
    if not _is_positive_finding(pf):
        return []
    names = []
    for part in pf.split("+"):
        part = part.strip()
        if not part:
            continue
        name = part.split("(")[0].split("[")[0].strip()
        if name:
            names.append(name)
    return names


def query_jobs_filtered(qs):
    q = (qs.get("q", [""])[0] or "").strip().lower()
    date_from = (qs.get("date_from", [""])[0] or "").strip()
    date_to = (qs.get("date_to", [""])[0] or "").strip()
    status = (qs.get("status", ["completed"])[0] or "completed").strip()
    specimen = (qs.get("specimen", ["all"])[0] or "all").strip()
    finding = (qs.get("finding", ["all"])[0] or "all").strip()
    page = max(1, int(qs.get("page", ["1"])[0] or 1))
    page_size = min(200, max(5, int(qs.get("page_size", ["25"])[0] or 25)))

    with DB_LOCK:
        conn = sqlite3.connect(JOBS_DB)
        conn.row_factory = dict_factory
        all_rows = conn.execute("SELECT * FROM jobs ORDER BY id DESC").fetchall()
        conn.close()

    specimens_avail = sorted({(r.get("specimen_type") or "unknown") for r in all_rows})
    matched = []
    for r in all_rows:
        if status != "all" and r.get("status") != status:
            continue
        if specimen != "all" and (r.get("specimen_type") or "unknown") != specimen:
            continue
        ts_ref = (r.get("finished_at") or r.get("created_at") or "")[:10]
        if date_from and ts_ref and ts_ref < date_from:
            continue
        if date_to and ts_ref and ts_ref > date_to:
            continue
        pf = r.get("primary_finding") or ""
        pos = _is_positive_finding(pf)
        if finding == "positive" and not pos:
            continue
        if finding == "negative" and pos:
            continue
        if finding not in ("all", "positive", "negative") and finding.lower() not in pf.lower():
            continue
        if q:
            hay = f"{r.get('sample_id','')} {pf} {r.get('specimen_type','')} {r.get('error_msg','')}".lower()
            if q not in hay:
                continue
        matched.append(r)

    dedup = qs.get("dedup", ["1"])[0] != "0"
    raw_matched_count = len(matched)
    if dedup:
        seen_sids = set()
        dedup_matched = []
        for r in matched:
            sid_key = (r.get("sample_id") or "").strip()
            if sid_key and sid_key in seen_sids:
                continue
            if sid_key:
                seen_sids.add(sid_key)
            dedup_matched.append(r)
        matched = dedup_matched
    dedup_removed = raw_matched_count - len(matched)

    # Aggregate statistics across ALL matched rows
    pathogen_counts = {}
    period_maps = {"day": {}, "week": {}, "month": {}}
    specimen_counts = {}
    completed_cnt = failed_cnt = pos_cnt = neg_cnt = 0
    for r in matched:
        st = r.get("status")
        if st == "completed":
            completed_cnt += 1
        elif st == "failed":
            failed_cnt += 1
        pf = r.get("primary_finding") or ""
        pos = _is_positive_finding(pf)
        if st == "completed":
            if pos:
                pos_cnt += 1
            else:
                neg_cnt += 1
        for name in _extract_pathogen_names(pf):
            pathogen_counts[name] = pathogen_counts.get(name, 0) + 1
        sp = r.get("specimen_type") or "unknown"
        specimen_counts[sp] = specimen_counts.get(sp, 0) + 1
        day = (r.get("finished_at") or r.get("created_at") or "unknown")[:10]
        try:
            d_obj = dt.date.fromisoformat(day)
            iso = d_obj.isocalendar()
            w_key = f"{iso[0]}-W{iso[1]:02d}"
            w_start = (d_obj - dt.timedelta(days=d_obj.weekday())).isoformat()
            w_end = (d_obj + dt.timedelta(days=6 - d_obj.weekday())).isoformat()
            m_key = day[:7]
            m_start = f"{m_key}-01"
            next_m = (d_obj.replace(day=28) + dt.timedelta(days=4))
            m_end = (next_m - dt.timedelta(days=next_m.day)).isoformat()
            keys = [("day", day, day, day), ("week", w_key, w_start, w_end), ("month", m_key, m_start, m_end)]
        except ValueError:
            keys = [("day", day, "", ""), ("week", day, "", ""), ("month", day, "", "")]
        for gran, pkey, p_from, p_to in keys:
            rec = period_maps[gran].setdefault(
                pkey, {"period": pkey, "date": pkey, "date_from": p_from, "date_to": p_to, "total": 0, "positive": 0, "negative": 0}
            )
            rec["total"] += 1
            if st == "completed":
                rec["positive" if pos else "negative"] += 1

    by_period = {}
    for gran, pmap in period_maps.items():
        lst = []
        for pkey in sorted(pmap.keys(), reverse=True):
            rec = pmap[pkey]
            done_n = rec["positive"] + rec["negative"]
            rec["positive_rate_pct"] = round(100.0 * rec["positive"] / max(1, done_n), 1) if done_n else 0.0
            lst.append(rec)
        by_period[gran] = lst

    top_pathogens = [
        {"name": k, "count": v, "pct": round(100.0 * v / max(1, completed_cnt), 1)}
        for k, v in sorted(pathogen_counts.items(), key=lambda kv: (-kv[1], kv[0]))[:15]
    ]
    by_date = by_period["day"]
    by_specimen = [{"name": k, "count": v} for k, v in sorted(specimen_counts.items(), key=lambda kv: -kv[1])]

    total_matched = len(matched)
    total_pages = max(1, (total_matched + page_size - 1) // page_size)
    page = min(page, total_pages)
    start = (page - 1) * page_size
    items = matched[start:start + page_size]

    return {
        "total_matched": total_matched,
        "page": page,
        "page_size": page_size,
        "total_pages": total_pages,
        "specimens_available": specimens_avail,
        "items": items,
        "stats": {
            "total_count": total_matched,
            "raw_matched_count": raw_matched_count,
            "dedup": dedup,
            "dedup_removed": dedup_removed,
            "completed_count": completed_cnt,
            "failed_count": failed_cnt,
            "positive_count": pos_cnt,
            "negative_count": neg_cnt,
            "positive_rate_pct": round(100.0 * pos_cnt / max(1, completed_cnt), 1),
            "top_pathogens": top_pathogens,
            "by_date": by_date,
            "by_period": by_period,
            "by_specimen": by_specimen,
        },
        "matched_all": matched,
    }


_PIE_COLORS = ["#2563eb", "#16a34a", "#d97706", "#7c3aed", "#0891b2", "#dc2626", "#ea580c", "#64748b"]


def _svg_pie_chart(items, title, unit="筆"):
    """Render an offline SVG donut pie chart with right-side legend."""
    total = sum(v for _, v in items)
    if not items or total <= 0:
        return f"<div class='chart-box'><h3>{html.escape(title)}</h3><p style='color:#64748b;font-size:12.5px;'>尚無資料</p></div>"

    cx, cy, r_out, r_in = 105, 105, 86, 44
    paths = []
    legend = []
    angle = -math.pi / 2.0

    for idx, (label, val) in enumerate(items):
        col = _PIE_COLORS[idx % len(_PIE_COLORS)]
        frac = val / total
        pct = frac * 100.0
        tip = f"{label}: {val} {unit} ({pct:.1f}%)"
        if len(items) == 1 or frac >= 0.9999:
            paths.append(
                f'<circle cx="{cx}" cy="{cy}" r="{(r_out + r_in) / 2:.1f}" fill="none" '
                f'stroke="{col}" stroke-width="{r_out - r_in:.1f}"><title>{html.escape(tip)}</title></circle>'
            )
        else:
            sweep = frac * 2.0 * math.pi
            a0, a1 = angle, angle + sweep
            angle = a1
            large = 1 if sweep > math.pi else 0
            x0, y0 = cx + r_out * math.cos(a0), cy + r_out * math.sin(a0)
            x1, y1 = cx + r_out * math.cos(a1), cy + r_out * math.sin(a1)
            xi1, yi1 = cx + r_in * math.cos(a1), cy + r_in * math.sin(a1)
            xi0, yi0 = cx + r_in * math.cos(a0), cy + r_in * math.sin(a0)
            d = (
                f"M {x0:.2f},{y0:.2f} A {r_out},{r_out} 0 {large} 1 {x1:.2f},{y1:.2f} "
                f"L {xi1:.2f},{yi1:.2f} A {r_in},{r_in} 0 {large} 0 {xi0:.2f},{yi0:.2f} Z"
            )
            paths.append(f'<path d="{d}" fill="{col}" stroke="#fff" stroke-width="1.5"><title>{html.escape(tip)}</title></path>')
        legend.append(
            f"<div class='leg-row'><span class='leg-dot' style='background:{col}'></span>"
            f"<span class='leg-name'>{html.escape(label)}</span>"
            f"<span class='leg-val'>{val} ({pct:.1f}%)</span></div>"
        )

    svg = (
        f'<svg viewBox="0 0 210 210" style="width:190px;height:190px;flex-shrink:0;">'
        f'{"".join(paths)}'
        f'<text x="{cx}" y="{cy - 4}" text-anchor="middle" font-size="18" font-weight="700" fill="#0f172a">{total}</text>'
        f'<text x="{cx}" y="{cy + 14}" text-anchor="middle" font-size="11" fill="#64748b">總計 ({html.escape(unit)})</text>'
        f'</svg>'
    )
    return (
        f"<div class='chart-box'><h3>{html.escape(title)}</h3>"
        f"<div class='pie-wrap'>{svg}<div class='leg-list'>{''.join(legend)}</div></div></div>"
    )


def _svg_pathogen_time_chart(dates, top_names, date_sp_counts, date_totals):
    """Render pathogen detection frequency (%) over time as SVG line or bar chart."""
    if not dates or not top_names:
        return "<p style='color:#64748b;font-size:12.5px;'>尚無特定病原隨時間檢出趨勢資料</p>"

    w, h = 920, 260
    pad_l, pad_r, pad_t, pad_b = 52, 20, 22, 42
    pw, ph = w - pad_l - pad_r, h - pad_t - pad_b
    grid = []
    for tick in (0, 25, 50, 75, 100):
        y = pad_t + ph * (1.0 - tick / 100.0)
        grid.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{w - pad_r}" y2="{y:.1f}" stroke="#e2e8f0" stroke-dasharray="3,3"/>')
        grid.append(f'<text x="{pad_l - 8}" y="{y + 4:.1f}" text-anchor="end" font-size="10.5" fill="#64748b">{tick}%</text>')

    series_svg = []
    legend = []
    if len(dates) == 1:
        d = dates[0]
        denom = max(1, date_totals.get(d, 0))
        bw = min(64.0, pw / max(1, len(top_names)) * 0.56)
        step = pw / max(1, len(top_names))
        for idx, name in enumerate(top_names):
            col = _PIE_COLORS[idx % len(_PIE_COLORS)]
            cnt = date_sp_counts.get(d, {}).get(name, 0)
            pct = cnt / denom * 100.0
            bh = ph * (pct / 100.0)
            x = pad_l + idx * step + (step - bw) / 2.0
            y = pad_t + ph - bh
            tip = f"{d} | {name}: {cnt}/{denom} ({pct:.1f}%)"
            series_svg.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{max(2.0, bh):.1f}" rx="4" fill="{col}">'
                f'<title>{html.escape(tip)}</title></rect>'
                f'<text x="{x + bw / 2:.1f}" y="{max(pad_t + 10, y - 5):.1f}" text-anchor="middle" font-size="10.5" font-weight="600" fill="#0f172a">{pct:.1f}% ({cnt})</text>'
                f'<text x="{x + bw / 2:.1f}" y="{h - 14}" text-anchor="middle" font-size="10.5" fill="#334155">{html.escape(name[:18])}</text>'
            )
            legend.append(f"<span class='leg-chip'><i style='background:{col}'></i>{html.escape(name)}</span>")
    else:
        n_d = len(dates)
        for i, d in enumerate(dates):
            x = pad_l + (i / (n_d - 1)) * pw
            if n_d <= 12 or i % max(1, n_d // 8) == 0 or i == n_d - 1:
                grid.append(f'<text x="{x:.1f}" y="{h - 14}" text-anchor="middle" font-size="10.5" fill="#475569">{html.escape(d)}</text>')
        for idx, name in enumerate(top_names):
            col = _PIE_COLORS[idx % len(_PIE_COLORS)]
            pts = []
            dots = []
            for i, d in enumerate(dates):
                denom = max(1, date_totals.get(d, 0))
                cnt = date_sp_counts.get(d, {}).get(name, 0)
                pct = cnt / denom * 100.0
                x = pad_l + (i / (n_d - 1)) * pw
                y = pad_t + ph * (1.0 - pct / 100.0)
                pts.append(f"{x:.1f},{y:.1f}")
                tip = f"{d} | {name}: {cnt}/{denom} ({pct:.1f}%)"
                dots.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{col}" stroke="#fff" stroke-width="1.5"><title>{html.escape(tip)}</title></circle>')
            series_svg.append(f'<polyline fill="none" stroke="{col}" stroke-width="2.4" points="{" ".join(pts)}"/>' + "".join(dots))
            legend.append(f"<span class='leg-chip'><i style='background:{col}'></i>{html.escape(name)}</span>")

    return (
        f"<div style='display:flex;flex-wrap:wrap;gap:10px;margin-bottom:8px;'>{''.join(legend)}</div>"
        f'<svg viewBox="0 0 {w} {h}" style="width:100%;height:auto;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">'
        f'{"".join(grid)}{"".join(series_svg)}</svg>'
    )


def render_cohort_html_report(res, qs):
    st = res["stats"]
    now_str = time.strftime("%Y-%m-%d %H:%M:%S %z")
    filter_desc = (
        f"狀態={qs.get('status',['completed'])[0]} | "
        f"起日={qs.get('date_from',['不限'])[0] or '不限'} | "
        f"迄日={qs.get('date_to',['不限'])[0] or '不限'} | "
        f"檢體={qs.get('specimen',['all'])[0]} | "
        f"結果={qs.get('finding',['all'])[0]} | "
        f"關鍵字={qs.get('q',['無'])[0] or '無'} | "
        f"去重={'僅採計最新 1 筆 (已排除 ' + str(st.get('dedup_removed', 0)) + ' 筆重複)' if st.get('dedup', True) else '否 (含重複紀錄)'}"
    )

    # Per-date pathogen frequency across matched rows
    date_sp_counts = {}
    date_totals = {}
    for r in res["matched_all"]:
        day = (r.get("finished_at") or r.get("created_at") or "unknown")[:10]
        if r.get("status") == "completed":
            date_totals[day] = date_totals.get(day, 0) + 1
            for name in _extract_pathogen_names(r.get("primary_finding") or ""):
                date_sp_counts.setdefault(day, {})[name] = date_sp_counts.get(day, {}).get(name, 0) + 1
    sorted_dates = sorted(date_totals.keys())
    top_names = [p["name"] for p in st["top_pathogens"][:6]]

    pie_pathogens = _svg_pie_chart([(p["name"], p["count"]) for p in st["top_pathogens"][:7]], "主要檢出病原體分佈 (Pie Chart)", "次")
    pie_specimens = _svg_pie_chart([(s["name"], s["count"]) for s in st.get("by_specimen", [])], "檢體類型分佈 (Pie Chart)", "筆")
    time_chart = _svg_pathogen_time_chart(sorted_dates, top_names, date_sp_counts, date_totals)

    freq_th = "".join(f"<th>{html.escape(n)}</th>" for n in top_names)
    freq_rows = []
    for d in reversed(sorted_dates):
        denom = max(1, date_totals.get(d, 0))
        cells = "".join(
            f"<td>{date_sp_counts.get(d, {}).get(n, 0)} ({date_sp_counts.get(d, {}).get(n, 0) / denom * 100:.1f}%)</td>"
            for n in top_names
        )
        freq_rows.append(f"<tr><td>{html.escape(d)}</td><td>{denom}</td>{cells}</tr>")
    freq_table = (
        f"<table><thead><tr><th>日期</th><th>已完成數</th>{freq_th}</tr></thead>"
        f"<tbody>{''.join(freq_rows) or '<tr><td colspan=\"2\">無資料</td></tr>'}</tbody></table>"
        if top_names else ""
    )

    path_rows = "".join(
        f"<tr><td>{html.escape(p['name'])}</td><td>{p['count']}</td><td>{p['pct']}%</td></tr>"
        for p in st["top_pathogens"]
    ) or "<tr><td colspan='3'>無檢出紀錄</td></tr>"
    date_rows = "".join(
        f"<tr><td>{html.escape(d['date'])}</td><td>{d['total']}</td><td>{d['positive']}</td><td>{d['negative']}</td></tr>"
        for d in st["by_date"]
    ) or "<tr><td colspan='4'>無資料</td></tr>"
    sample_rows = "".join(
        f"<tr><td>{r['id']}</td><td>{html.escape(r['sample_id'])}</td><td>{html.escape(r.get('specimen_type',''))}</td>"
        f"<td>{html.escape(r.get('status',''))}</td><td>{html.escape(r.get('created_at',''))}</td>"
        f"<td>{html.escape(r.get('finished_at',''))}</td><td>{html.escape(r.get('primary_finding') or r.get('error_msg') or '-')}</td></tr>"
        for r in res["matched_all"]
    )
    return f"""<!DOCTYPE html>
<html lang="zh-Hant"><head><meta charset="utf-8"><title>ClinMeta 批次統計報表</title>
<style>body{{font-family:system-ui,sans-serif;margin:24px;color:#0f172a;background:#fff}}h1{{margin:0 0 6px}}h2{{margin:22px 0 10px;font-size:17px;border-left:4px solid #2563eb;padding-left:8px}}h3{{margin:0 0 10px;font-size:14.5px;color:#1e293b}}
.meta{{color:#475569;font-size:13px;margin-bottom:18px}}.grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-bottom:20px}}
.charts-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:16px;margin-bottom:20px}}
.chart-box{{border:1px solid #cbd5e1;border-radius:10px;padding:14px;background:#f8fafc}}
.pie-wrap{{display:flex;align-items:center;gap:16px;flex-wrap:wrap}}
.leg-list{{flex:1;min-width:160px;font-size:12.5px}}.leg-row{{display:flex;align-items:center;gap:7px;margin:5px 0}}.leg-dot{{width:10px;height:10px;border-radius:999px;flex-shrink:0}}.leg-name{{flex:1;color:#1e293b}}.leg-val{{font-weight:600;color:#0f172a}}
.leg-chip{{display:inline-flex;align-items:center;gap:5px;font-size:12px;color:#334155}}.leg-chip i{{width:10px;height:10px;border-radius:2px;display:inline-block}}
.card{{border:1px solid #cbd5e1;border-radius:8px;padding:12px}}.card .v{{font-size:22px;font-weight:700;margin-top:4px}}
table{{border-collapse:collapse;width:100%;font-size:12.5px;margin-bottom:20px}}th,td{{border:1px solid #cbd5e1;padding:6px 8px;text-align:left}}
th{{background:#f1f5f9;cursor:pointer;user-select:none}}th:hover{{background:#e2e8f0}}</style></head><body>
<h1>ClinMeta 臨床宏基因組批次統計報表</h1>
<div class="meta">產出時間：{html.escape(now_str)}<br>篩選條件：{html.escape(filter_desc)}</div>
<div class="grid">
<div class="card"><div>符合筆數</div><div class="v">{st['total_count']}</div></div>
<div class="card"><div>已完成</div><div class="v">{st['completed_count']}</div></div>
<div class="card"><div>檢出陽性病原 (率)</div><div class="v">{st['positive_count']} ({st['positive_rate_pct']}%)</div></div>
<div class="card"><div>未檢出主要病原</div><div class="v">{st['negative_count']}</div></div>
</div>
<div class="charts-grid">{pie_pathogens}{pie_specimens}</div>
<h2>特定病原體隨時間檢出頻度趨勢（佔當日已完成樣本比例）</h2>
{time_chart}
{freq_table}
<h2>主要檢出病原統計 (Top Pathogens)</h2>
<table><thead><tr><th>病原物種</th><th>檢出樣本數</th><th>佔已完成比例</th></tr></thead><tbody>{path_rows}</tbody></table>
<h2>每日分析與陽性統計</h2>
<table><thead><tr><th>日期</th><th>總筆數</th><th>陽性數</th><th>未檢出數</th></tr></thead><tbody>{date_rows}</tbody></table>
<h2>符合條件樣本明細（共 {st['total_count']} 筆）</h2>
<table><thead><tr><th>ID</th><th>Sample ID</th><th>檢體</th><th>狀態</th><th>送件時間</th><th>完成時間</th><th>主要檢出摘要</th></tr></thead><tbody>{sample_rows}</tbody></table>
<script>
document.querySelectorAll('table').forEach(tbl=>{{
  const ths=tbl.querySelectorAll('thead th');
  ths.forEach((th,idx)=>{{
    th.dataset.orig=th.textContent.trim();
    th.onclick=()=>{{
      const tb=tbl.querySelector('tbody');if(!tb)return;
      const asc=th.dataset.dir!=='asc';
      ths.forEach(h=>{{h.dataset.dir='';if(h.dataset.orig)h.textContent=h.dataset.orig;}});
      th.dataset.dir=asc?'asc':'desc';th.textContent=th.dataset.orig+(asc?' ▲':' ▼');
      const pv=td=>{{const s=(td?td.textContent:'').trim(),c=s.split(' ')[0].replace(/[,%]/g,'');return /^-?\\d+(\\.\\d+)?$/.test(c)&&c!==''?{{n:parseFloat(c),ok:true,s}}:{{n:0,ok:false,s}};}};
      const rows=Array.from(tb.querySelectorAll('tr'));
      rows.sort((a,b)=>{{const va=pv(a.children[idx]),vb=pv(b.children[idx]);return(va.ok&&vb.ok)?(asc?va.n-vb.n:vb.n-va.n):(asc?va.s.localeCompare(vb.s,'zh-Hant'):vb.s.localeCompare(va.s,'zh-Hant'));}});
      rows.forEach(r=>tb.appendChild(r));
    }};
  }});
}});
</script></body></html>"""


def get_sample_result_bundle(sample_or_job_id):
    with DB_LOCK:
        conn = sqlite3.connect(JOBS_DB)
        conn.row_factory = dict_factory
        if str(sample_or_job_id).isdigit():
            job = conn.execute("SELECT * FROM jobs WHERE id = ?", (int(sample_or_job_id),)).fetchone()
        else:
            job = conn.execute("SELECT * FROM jobs WHERE sample_id = ? ORDER BY id DESC LIMIT 1", (sample_or_job_id,)).fetchone()
        conn.close()
    if not job:
        return None
    sample_id = job["sample_id"]

    outdir = job["outdir"]
    k2_html_fp = os.path.join(ROOT_DIR, "clean_results", "kraken2", "report.html") if sample_id == "sample_001" else os.path.join(outdir, "kraken2", "report.html")
    if not os.path.exists(k2_html_fp) and os.path.isdir(os.path.dirname(k2_html_fp)):
        k2_dir = os.path.dirname(k2_html_fp)
        if glob.glob(os.path.join(k2_dir, "*.report.txt")):
            k2_script = os.path.join(ROOT_DIR, "scripts", "kraken2_report_html.py")
            if os.path.exists(k2_script):
                subprocess.run([sys.executable, k2_script, k2_dir, k2_html_fp], check=False)
    files = {}
    candidates = {
        "sankey_html": os.path.join(outdir, "read_disposition_sankey.html"),
        "validation_html": os.path.join(outdir, "pathogen_validation_report.html"),
        "kraken2_html": k2_html_fp,
        "consensus_fasta": os.path.join(outdir, "influenza_a_h3n2_consensus.fasta"),
        "disposition_tsv": os.path.join(outdir, "read_disposition_summary.tsv"),
        "species_tsv": os.path.join(outdir, "species_validation_summary.tsv"),
        "unknown_tsv": os.path.join(outdir, "unknown_viral_screen.tsv"),
    }
    for k, fp in candidates.items():
        if os.path.exists(fp):
            rel = os.path.relpath(fp, ROOT_DIR)
            files[k] = "/files/" + urllib.parse.quote(rel)

    return {
        "job": job,
        "files": files,
        "disposition_rows": read_tsv(candidates["disposition_tsv"], 50),
        "species_rows": read_tsv(candidates["species_tsv"], 50),
        "unknown_rows": read_tsv(candidates["unknown_tsv"], 20),
    }


def _url_to_local_path(url_path):
    if not url_path or not url_path.startswith("/files/"):
        return None
    rel = urllib.parse.unquote(url_path[len("/files/"):]).lstrip("/")
    return os.path.abspath(os.path.join(ROOT_DIR, rel))


def _render_dict_table(rows, empty_msg="無資料"):
    if not rows:
        return f"<p class='muted'>{html.escape(empty_msg)}</p>"
    cols = list(rows[0].keys())
    th = "".join(f"<th>{html.escape(str(c))}</th>" for c in cols)
    trs = []
    for r in rows:
        tds = "".join(f"<td>{html.escape(str(r.get(c, '')))}</td>" for c in cols)
        trs.append(f"<tr>{tds}</tr>")
    return f"<table><thead><tr>{th}</tr></thead><tbody>{''.join(trs)}</tbody></table>"


def _parse_kraken2_top_species(outdir, limit=15):
    k2_dirs = []
    if outdir:
        k2_dirs.append(os.path.join(outdir, "kraken2"))
    k2_dirs.append(os.path.join(ROOT_DIR, "clean_results", "kraken2"))
    sections = []
    for k2_dir in k2_dirs:
        if not os.path.isdir(k2_dir):
            continue
        reports = sorted(glob.glob(os.path.join(k2_dir, "*.report*")))
        if not reports:
            continue
        for rp in reports[:4]:
            rows = []
            try:
                with open(rp, encoding="utf-8", errors="replace") as f:
                    for line in f:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 6:
                            continue
                        rank = parts[3].strip()
                        if rank != "S":
                            continue
                        pct = parts[0].strip()
                        clade = int(parts[1].strip())
                        direct = int(parts[2].strip())
                        taxid = parts[4].strip()
                        name = parts[5].strip()
                        if clade > 0:
                            rows.append({
                                "物種": name,
                                "TaxID": taxid,
                                "比例 (%)": pct,
                                "Clade Reads": f"{clade:,}",
                                "Direct Reads": f"{direct:,}",
                            })
            except Exception:
                continue
            rows.sort(key=lambda r: int(r["Clade Reads"].replace(",", "")), reverse=True)
            if rows:
                sections.append(
                    f"<h3>{html.escape(os.path.basename(rp))}（Top {min(limit, len(rows))} 物種）</h3>"
                    + _render_dict_table(rows[:limit])
                )
        if sections:
            break
    return "".join(sections) if sections else "<p class='muted'>尚無 Kraken2 物種分類報告明細。</p>"


def build_full_report_html(sample_or_job_id, auto_print=False):
    bundle = get_sample_result_bundle(sample_or_job_id)
    if not bundle:
        return None
    job = bundle["job"]
    files = bundle.get("files", {})
    sid = job.get("sample_id", "unknown")

    logo_b64 = ""
    logo_path = os.path.join(os.path.dirname(__file__), "assets", "clinmeta-logo-light.jpg")
    if os.path.exists(logo_path):
        with open(logo_path, "rb") as lf:
            logo_b64 = base64.b64encode(lf.read()).decode("ascii")
    logo_html = f'<img src="data:image/jpeg;base64,{logo_b64}" alt="ClinMeta" style="height:48px;border-radius:6px;">' if logo_b64 else ""

    sankey_svg = "<p class='muted'>無 Sankey 圖表資料。</p>"
    sankey_file = _url_to_local_path(files.get("sankey_html"))
    if sankey_file and os.path.isfile(sankey_file):
        with open(sankey_file, encoding="utf-8", errors="replace") as f:
            m = re.search(r"<svg\b.*?</svg>", f.read(), re.S)
            if m:
                sankey_svg = m.group(0)

    k2_html = _parse_kraken2_top_species(job.get("outdir"))
    auto_js = "window.addEventListener('load', () => setTimeout(() => window.print(), 350));" if auto_print else ""

    return f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<title>{html.escape(sid)}_ClinMeta_Full_Report</title>
<style>
  @page {{ size: A4; margin: 12mm; }}
  * {{ box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Noto Sans TC", sans-serif; margin: 0; background: #f8fafc; color: #0f172a; font-size: 12.5px; line-height: 1.45; -webkit-print-color-adjust: exact; print-color-adjust: exact; }}
  .toolbar {{ position: sticky; top: 0; z-index: 20; background: #0f172a; color: #fff; padding: 10px 20px; display: flex; align-items: center; justify-content: space-between; gap: 12px; }}
  .toolbar button {{ background: #2563eb; color: #fff; border: none; border-radius: 6px; padding: 7px 14px; font-weight: 600; cursor: pointer; font-size: 13px; }}
  .wrap {{ max-width: 1060px; margin: 18px auto; padding: 0 18px 28px; }}
  .card {{ background: #fff; border: 1px solid #e2e8f0; border-radius: 10px; padding: 16px 18px; margin-bottom: 14px; }}
  .hdr {{ display: flex; justify-content: space-between; align-items: center; gap: 16px; border-bottom: 2px solid #0f172a; padding-bottom: 12px; margin-bottom: 14px; }}
  h1 {{ font-size: 20px; margin: 0 0 4px; }}
  h2 {{ font-size: 15.5px; margin: 0 0 10px; border-left: 4px solid #2563eb; padding-left: 8px; color: #0f172a; }}
  h3 {{ font-size: 13.5px; margin: 12px 0 6px; color: #334155; }}
  .meta-grid {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; margin-bottom: 10px; }}
  .meta-item {{ background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 8px 10px; }}
  .meta-k {{ font-size: 11px; color: #64748b; font-weight: 600; }}
  .meta-v {{ font-size: 13px; font-weight: 700; margin-top: 2px; word-break: break-word; }}
  .finding-box {{ background: #f0fdf4; border: 1px solid #bbf7d0; border-left: 4px solid #16a34a; border-radius: 6px; padding: 10px 12px; font-weight: 600; color: #14532d; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 11.5px; margin-top: 6px; }}
  th, td {{ padding: 6px 8px; border-bottom: 1px solid #e2e8f0; text-align: left; }}
  th {{ background: #f1f5f9; font-weight: 600; color: #334155; }}
  .muted {{ color: #64748b; font-size: 12px; }}
  .page-break {{ page-break-before: always; }}
  @media print {{
    body {{ background: #fff; }}
    .no-print {{ display: none !important; }}
    .wrap {{ max-width: none; margin: 0; padding: 0; }}
    .card {{ border: 1px solid #cbd5e1; break-inside: avoid; }}
  }}
</style>
</head>
<body>
<div class="toolbar no-print">
  <span>ClinMeta 完整分析報告預覽（請於列印視窗選擇「<strong>另存為 PDF / Save as PDF</strong>」）</span>
  <button onclick="window.print()">列印 / 下載 PDF</button>
</div>
<div class="wrap">
  <div class="card">
    <div class="hdr">
      <div>
        <h1>ClinMeta 臨床微生物宏基因體完整分析報告</h1>
        <div class="muted">離線全流程（Hostile 去宿主 + Kraken2 PlusPF + 競爭性回貼驗證 + DIAMOND 病毒篩檢）</div>
      </div>
      {logo_html}
    </div>
    <div class="meta-grid">
      <div class="meta-item"><div class="meta-k">樣本編號 (Sample ID)</div><div class="meta-v">{html.escape(sid)}</div></div>
      <div class="meta-item"><div class="meta-k">檢體類型</div><div class="meta-v">{html.escape(job.get('specimen_type') or 'unknown')}</div></div>
      <div class="meta-item"><div class="meta-k">送件時間</div><div class="meta-v">{html.escape((job.get('created_at') or '').replace('T', ' ')[:19])}</div></div>
      <div class="meta-item"><div class="meta-k">完成時間</div><div class="meta-v">{html.escape((job.get('finished_at') or '').replace('T', ' ')[:19])}</div></div>
    </div>
    <div class="finding-box">主要檢出結論：{html.escape(job.get('primary_finding') or '尚無摘要')}</div>
  </div>

  <div class="card">
    <h2>一、全序列流向 Sankey 圖與互斥判定明細</h2>
    {sankey_svg}
    <h3>互斥歸屬明細表</h3>
    {_render_dict_table(bundle.get("disposition_rows", []), "無序列判定明細")}
  </div>

  <div class="card page-break">
    <h2>二、病原體競爭性回貼驗證總表</h2>
    {_render_dict_table(bundle.get("species_rows", []), "無顯著匹配病原體")}
  </div>

  <div class="card">
    <h2>三、Kraken2 PlusPF 物種分類摘要</h2>
    {k2_html}
  </div>

  <div class="card">
    <h2>四、未知／遠緣病毒篩檢（DIAMOND blastx）</h2>
    {_render_dict_table(bundle.get("unknown_rows", []), "無顯著遠緣／未知病毒訊號")}
  </div>
</div>
<script>{auto_js}</script>
</body>
</html>"""


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass

    def _send_json(self, obj, status=200):
        body = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_json(self):
        ln = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(ln) if ln > 0 else b"{}"
        return json.loads(raw.decode("utf-8"))

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path == "/" or path == "/index.html":
            idx = os.path.join(os.path.dirname(__file__), "index.html")
            with open(idx, "rb") as f:
                data = f.read()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return

        if path.startswith("/assets/"):
            rel_clean = posixpath.normpath(urllib.parse.unquote(path[len("/assets/"):])).lstrip("/")
            full = os.path.abspath(os.path.join(os.path.dirname(__file__), "assets", rel_clean))
            assets_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "assets"))
            if not full.startswith(assets_root) or not os.path.isfile(full):
                self.send_error(404)
                return
            ctype, _ = mimetypes.guess_type(full)
            with open(full, "rb") as f:
                data = f.read()
            self.send_response(200)
            self.send_header("Content-Type", ctype or "application/octet-stream")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return

        if path == "/api/health":
            self._send_json(check_offline_env())
            return

        if path == "/api/resources":
            with PROC_LOCK:
                active_count = len(ACTIVE_PROCS)
            self._send_json(compute_execution_plan(queued_count=1, active_count=active_count))
            return

        if path == "/api/pathogen-panel":
            cat_tsv = os.path.join(REF_DIR, "pathogen_panel", "PANEL_CATALOG.tsv")
            rows = read_tsv(cat_tsv, limit=500)
            self._send_json({"count": len(rows), "items": rows})
            return

        if path == "/api/input-scan":
            qs = urllib.parse.parse_qs(parsed.query)
            custom_dir = qs.get("dir", [""])[0]
            used_dir, samples = scan_input_samples(custom_dir)
            self._send_json({"input_dir": used_dir, "samples": samples})
            return

        if path == "/api/browse":
            qs = urllib.parse.parse_qs(parsed.query)
            req_path = qs.get("path", [ROOT_DIR])[0]
            self._send_json(browse_directory(req_path))
            return

        if path == "/api/jobs":
            self._send_json({"jobs": get_jobs()})
            return

        if path == "/api/results-query":
            qs = urllib.parse.parse_qs(parsed.query)
            res = query_jobs_filtered(qs)
            res.pop("matched_all", None)
            self._send_json(res)
            return

        if path == "/api/results-export":
            qs = urllib.parse.parse_qs(parsed.query)
            fmt = (qs.get("format", ["csv"])[0] or "csv").lower()
            res = query_jobs_filtered(qs)
            ts_tag = time.strftime("%Y%m%d_%H%M%S")
            if fmt == "html":
                inline = qs.get("inline", ["0"])[0] == "1"
                body = render_cohort_html_report(res, qs).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                if not inline:
                    self.send_header("Content-Disposition", f'attachment; filename="cohort_report_{ts_tag}.html"')
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            delim = "\t" if fmt == "tsv" else ","
            ext = "tsv" if fmt == "tsv" else "csv"
            buf = io.StringIO()
            w = csv.writer(buf, delimiter=delim)
            w.writerow(["id", "sample_id", "specimen_type", "status", "created_at", "started_at", "finished_at", "threads", "primary_finding", "error_msg", "outdir"])
            for r in res["matched_all"]:
                w.writerow([r["id"], r["sample_id"], r.get("specimen_type", ""), r.get("status", ""), r.get("created_at", ""), r.get("started_at", ""), r.get("finished_at", ""), r.get("threads", ""), r.get("primary_finding", ""), r.get("error_msg", ""), r.get("outdir", "")])
            body = ("\ufeff" + buf.getvalue()).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", f"text/{ext}; charset=utf-8")
            self.send_header("Content-Disposition", f'attachment; filename="cohort_report_{ts_tag}.{ext}"')
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        if path.startswith("/api/results/") and path.endswith("/full-report"):
            sid = urllib.parse.unquote(path[len("/api/results/"):-len("/full-report")])
            qs = urllib.parse.parse_qs(parsed.query)
            auto_print = qs.get("print", ["0"])[0] == "1"
            doc = build_full_report_html(sid, auto_print=auto_print)
            if not doc:
                self.send_error(404)
                return
            data = doc.encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return

        if path.startswith("/api/results/"):
            sid = urllib.parse.unquote(path[len("/api/results/"):])
            bundle = get_sample_result_bundle(sid)
            if not bundle:
                self._send_json({"error": "Sample not found"}, 404)
            else:
                self._send_json(bundle)
            return

        if path.startswith("/files/"):
            rel = urllib.parse.unquote(path[len("/files/"):])
            rel_clean = posixpath.normpath(rel).lstrip("/")
            full = os.path.abspath(os.path.join(ROOT_DIR, rel_clean))
            allowed = (
                full.startswith(os.path.abspath(OUTPUT_DIR))
                or full.startswith(os.path.abspath(os.path.join(ROOT_DIR, "analysis")))
                or full.startswith(os.path.abspath(os.path.join(ROOT_DIR, "clean_results")))
                or full.startswith(os.path.abspath(os.path.join(REF_DIR, "pathogen_panel")))
            )
            if not allowed or not os.path.isfile(full):
                self.send_error(404)
                return
            ctype, _ = mimetypes.guess_type(full)
            with open(full, "rb") as f:
                data = f.read()
            self.send_response(200)
            self.send_header("Content-Type", (ctype or "application/octet-stream") + ("; charset=utf-8" if ctype and ctype.startswith("text/") else ""))
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return

        self.send_error(404)

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path == "/api/jobs":
            data = self._read_json()
            items = data.get("samples") or [data]
            created = []
            now = time.strftime("%Y-%m-%dT%H:%M:%S%z")
            ts_tag = time.strftime("%Y%m%d_%H%M%S")
            with DB_LOCK:
                conn = sqlite3.connect(JOBS_DB)
                max_prio = conn.execute("SELECT COALESCE(MAX(priority), 10) FROM jobs").fetchone()[0]
                for idx, it in enumerate(items, 1):
                    sid = (it.get("sample_id") or "").strip()
                    if not sid:
                        continue
                    suffix = f"{ts_tag}_{idx:02d}" if len(items) > 1 else ts_tag
                    outdir = os.path.join(OUTPUT_DIR, f"{sid}_{suffix}")
                    cur = conn.execute(
                        """
                        INSERT INTO jobs (sample_id, dna_path, cdna_path, specimen_type, threads, priority, status, stage, outdir, created_at)
                        VALUES (?, ?, ?, ?, ?, ?, 'queued', 'Queued', ?, ?)
                        """,
                        (
                            sid,
                            it.get("dna_path", ""),
                            it.get("cdna_path", ""),
                            it.get("specimen_type", "unknown"),
                            int(it.get("threads", 16)),
                            max_prio + idx,
                            outdir,
                            now,
                        ),
                    )
                    created.append(cur.lastrowid)
                conn.commit()
                conn.close()
            self._send_json({"queued_ids": created}, 201)
            return

        if path == "/api/jobs/reorder":
            data = self._read_json()
            jid = int(data["id"])
            direction = data.get("direction", "up")
            with DB_LOCK:
                conn = sqlite3.connect(JOBS_DB)
                conn.row_factory = dict_factory
                queued = conn.execute("SELECT id, priority FROM jobs WHERE status='queued' ORDER BY priority ASC, id ASC").fetchall()
                idx = next((i for i, x in enumerate(queued) if x["id"] == jid), None)
                if idx is not None:
                    swap_idx = idx - 1 if direction == "up" else idx + 1
                    if 0 <= swap_idx < len(queued):
                        a, b = queued[idx], queued[swap_idx]
                        conn.execute("UPDATE jobs SET priority=? WHERE id=?", (b["priority"], a["id"]))
                        conn.execute("UPDATE jobs SET priority=? WHERE id=?", (a["priority"], b["id"]))
                        conn.commit()
                conn.close()
            self._send_json({"ok": True})
            return

        if path.startswith("/api/jobs/") and path.endswith("/cancel"):
            jid = int(path.split("/")[3])
            with DB_LOCK:
                conn = sqlite3.connect(JOBS_DB)
                conn.execute("UPDATE jobs SET status='cancelled', stage='Cancelled' WHERE id=? AND status='queued'", (jid,))
                conn.commit()
                conn.close()
            self._send_json({"ok": True})
            return

        if path.startswith("/api/jobs/") and path.endswith("/retry"):
            jid = int(path.split("/")[3])
            with DB_LOCK:
                conn = sqlite3.connect(JOBS_DB)
                conn.execute(
                    "UPDATE jobs SET status='queued', stage='Queued', progress_pct=0, error_msg='', finished_at='' WHERE id=? AND status IN ('failed', 'cancelled')",
                    (jid,),
                )
                conn.commit()
                conn.close()
            self._send_json({"ok": True})
            return

        if path.startswith("/api/jobs/") and path.endswith("/delete"):
            jid = int(path.split("/")[3])
            with DB_LOCK:
                conn = sqlite3.connect(JOBS_DB)
                conn.row_factory = dict_factory
                job = conn.execute("SELECT * FROM jobs WHERE id=?", (jid,)).fetchone()
                if not job:
                    conn.close()
                    self._send_json({"error": "Job not found"}, 404)
                    return
                if job.get("status") == "running":
                    conn.close()
                    self._send_json({"error": "Cannot delete a running job"}, 400)
                    return
                outdir = os.path.abspath(job.get("outdir") or "")
                conn.execute("DELETE FROM jobs WHERE id=?", (jid,))
                conn.commit()
                other_ref = conn.execute("SELECT COUNT(*) AS c FROM jobs WHERE outdir=?", (job.get("outdir"),)).fetchone()["c"]
                conn.close()
            out_root = os.path.abspath(OUTPUT_DIR)
            if other_ref == 0 and outdir.startswith(out_root + os.sep) and os.path.isdir(outdir):
                shutil.rmtree(outdir, ignore_errors=True)
            self._send_json({"ok": True, "deleted_id": jid})
            return

        self.send_error(404)


def main():
    init_db()
    warm_kraken2_cache(DB_DIR, async_run=True)
    t = threading.Thread(target=worker_loop, daemon=True)
    t.start()
    srv = ThreadingHTTPServer(("0.0.0.0", PORT), Handler)
    print(f"CDC Offline Pathogen Web GUI listening on http://0.0.0.0:{PORT}", flush=True)
    srv.serve_forever()


if __name__ == "__main__":
    main()

