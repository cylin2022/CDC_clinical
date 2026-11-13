#!/usr/bin/env python3
"""Fully offline single-sample Nanopore pathogen analysis runner.

Takes DNA and/or cDNA FASTQ inputs, runs:
  1) chopper QC + hostile dehosting (using local human-masked indices)
  2) Kraken2 PlusPF classification
  3) Offline competitive minimap2 validation against human-masked + matched pathogen_panel refs
  4) Unknown/divergent virus DIAMOND blastx screen on non-TSO unclassified reads
  5) Generates species_validation_summary.tsv, read_disposition_summary.tsv,
     read_disposition_sankey.html, and pathogen_validation_report.html.
"""

import argparse
import base64
import collections
from concurrent.futures import ThreadPoolExecutor
import csv
import glob
import gzip
import html
import json
import math
import math
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from resource_optimizer import compute_execution_plan  # noqa: E402

TSO_MOTIFS = (
    "GTATCAACGCAGAGT",
    "ACTCTGCGTTGATAC",
    "AGTGGTATCAAC",
    "GTTGATACCACT",
    "AAGCAGTGGTAT",
    "ATACCACTGCTT",
)


def update_progress(outdir, stage, pct, message, extra=None):
    os.makedirs(outdir, exist_ok=True)
    payload = {
        "stage": stage,
        "pct": pct,
        "message": message,
        "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    if extra:
        payload.update(extra)
    with open(os.path.join(outdir, "progress.json"), "w") as f:
        json.dump(payload, f, indent=2)


def find_fastqs(path_or_dir):
    if not path_or_dir:
        return []
    if os.path.isfile(path_or_dir):
        return [path_or_dir]
    if os.path.isdir(path_or_dir):
        return sorted(
            glob.glob(os.path.join(path_or_dir, "*.fastq.gz"))
            + glob.glob(os.path.join(path_or_dir, "*.fq.gz"))
            + glob.glob(os.path.join(path_or_dir, "*.fastq"))
        )
    return []


def count_fastq_reads(fq_path):
    n = 0
    opener = gzip.open if fq_path.endswith(".gz") else open
    with opener(fq_path, "rt") as f:
        for i, _ in enumerate(f, 1):
            if i % 4 == 0:
                n += 1
    return n


def run_qc_and_dehost(fastqs, lib_label, preset, index_prefix, out_fastq_gz, threads, tmp_dir):
    """Run chopper QC and hostile clean on FASTQs for one library."""
    os.makedirs(tmp_dir, exist_ok=True)
    merged_raw = os.path.join(tmp_dir, f"{lib_label}.raw.fastq.gz")
    if len(fastqs) == 1 and fastqs[0].endswith(".gz"):
        merged_raw = fastqs[0]
    else:
        with open(merged_raw, "wb") as out:
            for fq in fastqs:
                with open(fq, "rb") as inp:
                    while True:
                        buf = inp.read(1024 * 1024)
                        if not buf:
                            break
                        out.write(buf)

    raw_reads = count_fastq_reads(merged_raw)
    qc_fq = os.path.join(tmp_dir, f"{lib_label}.qc.fastq.gz")
    chop_cmd = (
        f"zcat -f {merged_raw!r} | chopper -q 10 -l 200 --headcrop 15 --tailcrop 15 "
        f"--threads {threads} | pigz -p 4 > {qc_fq!r}"
    )
    subprocess.run(["bash", "-c", chop_cmd], check=True)
    qc_reads = count_fastq_reads(qc_fq)

    hostile_dir = os.path.join(tmp_dir, f"hostile_{lib_label}")
    os.makedirs(hostile_dir, exist_ok=True)
    for old in glob.glob(os.path.join(hostile_dir, "*.fastq.gz")):
        os.remove(old)

    if qc_reads > 0:
        subprocess.run(
            [
                "hostile", "clean",
                "--fastq1", qc_fq,
                "--aligner", "minimap2",
                "--index", index_prefix,
                f"--aligner-args=-x {preset}",
                "--airplane",
                "--threads", str(threads),
                "--output", hostile_dir,
                "--force",
            ],
            check=True,
        )
        cleaned = glob.glob(os.path.join(hostile_dir, "*.fastq.gz"))
        if not cleaned:
            raise RuntimeError(f"hostile produced no output in {hostile_dir}")
        os.replace(cleaned[0], out_fastq_gz)
    else:
        with gzip.open(out_fastq_gz, "wb"):
            pass

    nonhost_reads = count_fastq_reads(out_fastq_gz)
    return {
        "raw_reads": raw_reads,
        "qc_reads": qc_reads,
        "qc_failed": raw_reads - qc_reads,
        "host_removed": qc_reads - nonhost_reads,
        "nonhost_reads": nonhost_reads,
    }


def select_panel_refs(k2_reports, panel_dir, max_refs=12):
    """Select top species taxids from Kraken2 reports that exist in offline pathogen_panel."""
    panel_by_taxid = {}
    for fp in sorted(glob.glob(os.path.join(panel_dir, "*__*.fna.gz"))):
        tid = os.path.basename(fp).split("__", 1)[0]
        panel_by_taxid.setdefault(tid, []).append(fp)

    scores = collections.Counter()
    names = {}
    for rep in k2_reports:
        if not os.path.exists(rep):
            continue
        with open(rep) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                clade_reads = int(parts[1])
                dmin = int(parts[4])
                rank = parts[5]
                tid = parts[6]
                name = parts[7].strip()
                if tid != "9606" and rank in ("S", "S1") and tid in panel_by_taxid and clade_reads >= 3:
                    scores[tid] += dmin + clade_reads
                    names[tid] = name

    selected = []
    for tid, _ in scores.most_common(max_refs):
        selected.append((tid, names[tid], panel_by_taxid[tid][0]))
    return selected


def run_competitive_validation(libs_clean, selected_refs, human_fa_gz, outdir, threads, use_shm_tmp=False):
    """Build sample-specific offline competitive reference (human + matched panel refs) and map."""
    if not selected_refs:
        return [], {}

    work_tmp = os.path.join("/dev/shm", f"cdc_comp_{os.getpid()}") if use_shm_tmp else os.path.join(outdir, "tmp")
    os.makedirs(work_tmp, exist_ok=True)
    comp_fa = os.path.join(work_tmp, "sample_competitive.fa")
    comp_mmi = os.path.join(work_tmp, "sample_competitive.mmi")
    if True:
        with open(comp_fa, "wt") as out:
            with gzip.open(human_fa_gz, "rt") as inp:
                for line in inp:
                    if line.startswith(">"):
                        cid = line[1:].strip().split()[0]
                        out.write(f">HUMAN|{cid}\n")
                    else:
                        out.write(line)
            for tid, name, fna in selected_refs:
                tag = f"TAX_{tid}"
                with gzip.open(fna, "rt") as inp:
                    for line in inp:
                        if line.startswith(">"):
                            cid = line[1:].strip().split()[0]
                            out.write(f">{tag}|{cid}\n")
                        else:
                            out.write(line)

        subprocess.run(["minimap2", "-x", "map-ont", "-t", str(threads), "-d", comp_mmi, comp_fa], check=True)

    ref_names = {f"TAX_{tid}": name for tid, name, _ in selected_refs}
    species_rows = []
    map_overview = {}

    for lib_label, fq in libs_clean:
        if not os.path.exists(fq):
            continue
        bam = os.path.join(outdir, f"{lib_label}.nonhuman.bam")
        unsorted = bam + ".unsorted.bam"
        mm2 = subprocess.Popen(
            ["minimap2", "-ax", "map-ont", "--secondary=no", "-t", str(threads), comp_mmi, fq],
            stdout=subprocess.PIPE, text=True,
        )
        sv = subprocess.Popen(
            ["samtools", "view", "-b", "-@", "4", "-o", unsorted, "-"],
            stdin=subprocess.PIPE, text=True,
        )
        counts = {"total": 0, "unmapped": 0, "human": 0, "nonhuman": 0}
        sp_reads = collections.Counter()
        sp_q20 = collections.Counter()

        for line in mm2.stdout:
            if line.startswith("@"):
                if line.startswith("@SQ") and "\tSN:HUMAN|" in line:
                    continue
                sv.stdin.write(line)
                continue
            fields = line.split("\t", 6)
            flag = int(fields[1])
            if flag & 0x900:
                continue
            counts["total"] += 1
            if flag & 0x4:
                counts["unmapped"] += 1
                continue
            rname = fields[2]
            if rname.startswith("HUMAN|"):
                counts["human"] += 1
                continue
            counts["nonhuman"] += 1
            tag = rname.split("|", 1)[0]
            sp_reads[tag] += 1
            if int(fields[4]) >= 20:
                sp_q20[tag] += 1
            sv.stdin.write(line)

        sv.stdin.close()
        mm2.wait()
        sv.wait()
        subprocess.run(["samtools", "sort", "-@", "4", "-o", bam, unsorted], check=True)
        subprocess.run(["samtools", "index", bam], check=True)
        os.remove(unsorted)

        cov_out = subprocess.check_output(["samtools", "coverage", bam], text=True)
        sp_cov = collections.defaultdict(lambda: {"ref_len": 0, "cov_bases": 0, "depth_sum": 0.0})
        for line in cov_out.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            tag = cols[0].split("|", 1)[0]
            ln = int(cols[2])
            cb = int(cols[4])
            md = float(cols[6])
            sp_cov[tag]["ref_len"] += ln
            sp_cov[tag]["cov_bases"] += cb
            sp_cov[tag]["depth_sum"] += md * ln

        map_overview[lib_label] = counts
        for tid, name, _ in selected_refs:
            tag = f"TAX_{tid}"
            cv = sp_cov[tag]
            rlen = cv["ref_len"]
            breadth = (cv["cov_bases"] / rlen * 100.0) if rlen else 0.0
            mdepth = (cv["depth_sum"] / rlen) if rlen else 0.0
            species_rows.append({
                "library": lib_label,
                "taxid": tid,
                "species": name,
                "ref_len": rlen,
                "mapped_reads": sp_reads[tag],
                "mapped_reads_mapq20": sp_q20[tag],
                "covered_bases": cv["cov_bases"],
                "breadth_pct": round(breadth, 4),
                "mean_depth": round(mdepth, 4),
            })

    if os.path.exists(comp_fa):
        os.remove(comp_fa)
    if os.path.exists(comp_mmi):
        os.remove(comp_mmi)
    if use_shm_tmp and os.path.isdir(work_tmp):
        shutil.rmtree(work_tmp, ignore_errors=True)

    tsv_path = os.path.join(outdir, "species_validation_summary.tsv")
    if species_rows:
        with open(tsv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(species_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(species_rows)
    return species_rows, map_overview


def _load_clinmeta_logo_html():
    """Return base64-embedded ClinMeta logo img tag if available locally."""
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    logo_path = os.path.join(repo_root, "webapp", "assets", "clinmeta-logo-light.jpg")
    if os.path.exists(logo_path):
        with open(logo_path, "rb") as lf:
            b64 = base64.b64encode(lf.read()).decode("ascii")
        return f'<img src="data:image/jpeg;base64,{b64}" alt="ClinMeta" style="height:52px;border-radius:6px;">'
    return ""


def _render_sankey_svg(nodes, links, total_all, width=1160, height=540, title=""):
    """Render a multi-column SVG Sankey diagram with power=0.45 scaling."""
    cols = sorted({n["col"] for n in nodes.values()})
    if not cols:
        return ""
    col_x = {c: 24 + i * ((width - 330) / max(1, len(cols) - 1)) for i, c in enumerate(cols)}
    node_w = 18
    pad_top = 46
    pad_bot = 22
    avail_h = height - pad_top - pad_bot

    for nid, n in nodes.items():
        in_v = sum(l["val"] for l in links if l["target"] == nid)
        out_v = sum(l["val"] for l in links if l["source"] == nid)
        n["val"] = max(in_v, out_v, n.get("val", 0))

    def scale_v(v):
        return math.pow(max(v, 1), 0.45)

    for c in cols:
        cnodes = [nid for nid, n in nodes.items() if n["col"] == c]
        gap = 12
        total_scaled = sum(scale_v(nodes[nid]["val"]) for nid in cnodes) or 1.0
        usable_h = max(40.0, avail_h - gap * max(0, len(cnodes) - 1))
        y = pad_top
        for nid in cnodes:
            nh = max(6.0, usable_h * (scale_v(nodes[nid]["val"]) / total_scaled))
            nodes[nid]["x"] = col_x[c]
            nodes[nid]["y"] = y
            nodes[nid]["h"] = nh
            nodes[nid]["out_y"] = y
            nodes[nid]["in_y"] = y
            y += nh + gap

    out_totals = {nid: sum(l["val"] for l in links if l["source"] == nid) for nid in nodes}
    in_totals = {nid: sum(l["val"] for l in links if l["target"] == nid) for nid in nodes}

    svg_links = []
    for l in links:
        s = nodes[l["source"]]
        t = nodes[l["target"]]
        sh = s["h"] * (l["val"] / max(1, out_totals[l["source"]]))
        th = t["h"] * (l["val"] / max(1, in_totals[l["target"]]))
        y0 = s["out_y"] + sh / 2.0
        y1 = t["in_y"] + th / 2.0
        s["out_y"] += sh
        t["in_y"] += th

        x0 = s["x"] + node_w
        x1 = t["x"]
        xm = (x0 + x1) / 2.0
        stroke_w = max(2.0, (sh + th) / 2.0)
        color = l.get("color", "#94a3b8")
        pct = l["val"] / max(1, total_all) * 100.0
        tip = f"{s['label']} -> {t['label']}: {l['val']:,} reads ({pct:.3f}%)"
        svg_links.append(
            f'<path d="M {x0:.1f},{y0:.1f} C {xm:.1f},{y0:.1f} {xm:.1f},{y1:.1f} {x1:.1f},{y1:.1f}" '
            f'fill="none" stroke="{color}" stroke-opacity="0.38" stroke-width="{stroke_w:.1f}" class="sankey-link">'
            f'<title>{html.escape(tip)}</title></path>'
        )

    svg_nodes = []
    for nid, n in nodes.items():
        x, y, h = n["x"], n["y"], n["h"]
        color = n.get("color", "#334155")
        lbl = f"{n['label']} ({n['val']:,})"
        tx = x + node_w + 6
        svg_nodes.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{node_w}" height="{h:.1f}" rx="3" fill="{color}">'
            f'<title>{html.escape(lbl)}</title></rect>'
            f'<text x="{tx:.1f}" y="{y + h / 2 + 4:.1f}" font-size="11.5" font-weight="600" fill="#0f172a" text-anchor="start">{html.escape(lbl)}</text>'
        )

    return f"""<svg viewBox="0 0 {width} {height}" style="width:100%;height:auto;background:#fff;border:1px solid #e2e8f0;border-radius:10px;">
  <text x="20" y="26" font-size="14" font-weight="700" fill="#0f172a">{html.escape(title)}</text>
  {''.join(svg_links)}
  {''.join(svg_nodes)}
</svg>"""


def write_simple_reports(sample_id, dehost_stats, species_rows, map_overview, outdir):
    """Write read_disposition_summary.tsv, read_disposition_sankey.html, and pathogen_validation_report.html."""
    disp_rows = []
    total_raw = sum(s["raw_reads"] for s in dehost_stats.values())
    lib_totals = {lib: max(1, st["raw_reads"]) for lib, st in dehost_stats.items()}
    for lib, st in dehost_stats.items():
        ov = map_overview.get(lib, {"human": 0, "nonhuman": 0, "unmapped": st["nonhost_reads"]})
        disp_rows.append((lib, "QC_Failed", "chopper", "Filtered out", st["qc_failed"]))
        disp_rows.append((lib, "Host_Removed", "hostile", "Host (Homo sapiens)", st["host_removed"]))
        disp_rows.append((lib, "NonHost_Residual_Human", "Competitive minimap2", "Residual Host", ov.get("human", 0)))
        disp_rows.append((lib, "NonHost_Pathogen_Mapped", "Competitive minimap2", "Mapped to Pathogen Panel", ov.get("nonhuman", 0)))
        disp_rows.append((lib, "NonHost_Unmapped", "Kraken2 + Offline Screen", "Background / Primer / Unclassified", ov.get("unmapped", 0)))

    disp_tsv = os.path.join(outdir, "read_disposition_summary.tsv")
    with open(disp_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["library", "final_category", "method", "interpretation", "reads", "pct_of_total_raw"])
        for lib, cat, m, interp, cnt in disp_rows:
            w.writerow([lib, cat, m, interp, cnt, f"{cnt / max(1, total_raw) * 100:.4f}"])

    logo_html = _load_clinmeta_logo_html()
    base_css = """
      body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; background: #f1f5f9; color: #0f172a; line-height: 1.5; }
      .container { max-width: 1200px; margin: 20px auto; padding: 0 18px; }
      .header-card { display: flex; align-items: center; justify-content: space-between; gap: 16px; flex-wrap: wrap; margin-bottom: 16px; background: #fff; padding: 14px 20px; border-radius: 10px; box-shadow: 0 1px 3px rgba(15,23,42,0.08); }
      .header-card h1 { font-size: 21px; margin: 0 0 4px; }
      .sub { color: #475569; font-size: 13.5px; margin: 0; }
      .panel { background: #fff; border-radius: 10px; padding: 18px 20px; margin-bottom: 18px; box-shadow: 0 1px 3px rgba(15,23,42,0.08); }
      .panel h2 { font-size: 16.5px; margin: 0 0 12px; border-bottom: 1px solid #e2e8f0; padding-bottom: 8px; color: #0f172a; }
      .cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 14px; margin-bottom: 18px; }
      .card { background: #fff; border-radius: 10px; padding: 14px 16px; box-shadow: 0 1px 3px rgba(15,23,42,0.08); border-left: 4px solid #2563eb; }
      .card.green { border-left-color: #16a34a; }
      .card.amber { border-left-color: #d97706; }
      .card.slate { border-left-color: #64748b; }
      .card-title { font-size: 11.5px; text-transform: uppercase; letter-spacing: 0.04em; color: #64748b; font-weight: 600; }
      .card-val { font-size: 20px; font-weight: 700; margin: 4px 0; color: #0f172a; }
      .card-note { font-size: 12.5px; color: #475569; }
      table { width: 100%; border-collapse: collapse; font-size: 13px; }
      th, td { padding: 8px 10px; border-bottom: 1px solid #e2e8f0; text-align: left; }
      th { background: #f8fafc; font-weight: 600; color: #334155; cursor: pointer; user-select: none; }
      th:hover { background: #e2e8f0; }
      td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
      code { background: #f1f5f9; padding: 2px 5px; border-radius: 4px; font-size: 12px; }
      .badge { display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11px; font-weight: 600; margin-left: 6px; }
      .badge-pos { background: #dcfce7; color: #166534; }
      .badge-warn { background: #fef3c7; color: #92400e; }
      .badge-neg { background: #e2e8f0; color: #334155; }
      .bar-row { display: flex; align-items: center; gap: 10px; margin: 8px 0; font-size: 13px; }
      .bar-label { width: 250px; color: #1e293b; font-weight: 500; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
      .bar-track { flex: 1; height: 14px; background: #e2e8f0; border-radius: 7px; overflow: hidden; }
      .bar-fill { height: 100%; background: linear-gradient(90deg, #16a34a, #22c55e); border-radius: 7px; }
      .bar-fill.warn { background: linear-gradient(90deg, #d97706, #f59e0b); }
      .bar-val { width: 220px; text-align: right; font-variant-numeric: tabular-nums; color: #334155; font-size: 12.5px; }
      .sankey-link:hover { stroke-opacity: 0.78; }
    """

    # --- 1. Tab 1: Interactive Multi-Column Sankey Diagram + Mutually-Exclusive Disposition Table ---
    nodes = {"raw": {"col": 0, "label": "All Raw Reads", "color": "#1e293b"}}
    links = []
    lib_colors = {"DNA": "#2563eb", "cDNA": "#7c3aed"}

    tot_qc = sum(st["qc_failed"] for st in dehost_stats.values())
    tot_host = sum(st["host_removed"] for st in dehost_stats.values())
    tot_nonhost = sum(st["nonhost_reads"] for st in dehost_stats.values())
    if tot_qc > 0:
        nodes["qc_fail"] = {"col": 2, "label": "QC Failed (chopper)", "color": "#94a3b8"}
    if tot_host > 0:
        nodes["host"] = {"col": 2, "label": "Host Removed (hostile)", "color": "#0284c7"}
    if tot_nonhost > 0:
        nodes["nonhost"] = {"col": 2, "label": "Non-Host Retained", "color": "#0f766e"}

    for lib, st in dehost_stats.items():
        lid = f"lib_{lib}"
        lcol = lib_colors.get(lib, "#3b82f6")
        nodes[lid] = {"col": 1, "label": f"{lib} Library", "color": lcol}
        if st["raw_reads"] > 0:
            links.append({"source": "raw", "target": lid, "val": st["raw_reads"], "color": lcol})
        if st["qc_failed"] > 0:
            links.append({"source": lid, "target": "qc_fail", "val": st["qc_failed"], "color": "#94a3b8"})
        if st["host_removed"] > 0:
            links.append({"source": lid, "target": "host", "val": st["host_removed"], "color": "#0284c7"})
        if st["nonhost_reads"] > 0:
            links.append({"source": lid, "target": "nonhost", "val": st["nonhost_reads"], "color": "#0f766e"})

    tot_res_human = sum(map_overview.get(lib, {}).get("human", 0) for lib in dehost_stats)
    tot_pathogen = sum(map_overview.get(lib, {}).get("nonhuman", 0) for lib in dehost_stats)
    tot_unmapped = sum(map_overview.get(lib, {}).get("unmapped", 0) for lib in dehost_stats)

    if tot_nonhost > 0:
        if tot_res_human > 0:
            nodes["res_human"] = {"col": 3, "label": "Residual Human (mapped)", "color": "#0ea5e9"}
            links.append({"source": "nonhost", "target": "res_human", "val": tot_res_human, "color": "#0ea5e9"})

        # Break down pathogen-mapped reads by species (proportional if sum differs slightly)
        sp_totals = collections.Counter()
        for r in species_rows:
            sp_totals[r["species"]] += int(r["mapped_reads"])
        assigned_sp = 0
        if tot_pathogen > 0 and sp_totals:
            top_sp = sp_totals.most_common(5)
            sum_top = sum(v for _, v in top_sp)
            for idx, (sp_name, sp_cnt) in enumerate(top_sp):
                if sp_cnt <= 0:
                    continue
                # Clamp so total pathogen links equal tot_pathogen
                val = min(sp_cnt, tot_pathogen - assigned_sp)
                if val <= 0:
                    continue
                assigned_sp += val
                nid = f"sp_{idx}"
                col_hex = "#16a34a" if idx == 0 else ("#15803d" if idx == 1 else "#0d9488")
                nodes[nid] = {"col": 3, "label": sp_name, "color": col_hex}
                links.append({"source": "nonhost", "target": nid, "val": val, "color": col_hex})
        rem_path = tot_pathogen - assigned_sp
        if rem_path > 0:
            nodes["sp_other"] = {"col": 3, "label": "Other Pathogen Panel", "color": "#f59e0b"}
            links.append({"source": "nonhost", "target": "sp_other", "val": rem_path, "color": "#f59e0b"})

        if tot_unmapped > 0:
            nodes["unmapped"] = {"col": 3, "label": "Unmapped / Background / Primer", "color": "#d97706"}
            links.append({"source": "nonhost", "target": "unmapped", "val": tot_unmapped, "color": "#d97706"})

    lib_desc = " + ".join(f"{html.escape(lib)} <code>{st['raw_reads']:,}</code>" for lib, st in dehost_stats.items())
    sankey_svg = _render_sankey_svg(
        nodes,
        links,
        total_raw,
        width=1160,
        height=520,
        title=f"{sample_id} 全序列流向與最終判定 Sankey 圖（總計 {total_raw:,} Raw Reads；採冪次縮放以兼顧微量病原體可視性）",
    )

    d_html_rows = []
    for lib, cat, m, interp, cnt in disp_rows:
        hl = " style='background:#f0fdf4;font-weight:600;'" if cat == "NonHost_Pathogen_Mapped" and cnt > 0 else ""
        pct_lib = cnt / lib_totals[lib] * 100.0
        pct_all = cnt / max(1, total_raw) * 100.0
        d_html_rows.append(
            f"<tr{hl}><td>{html.escape(lib)}</td><td><code>{html.escape(cat)}</code></td><td>{html.escape(m)}</td>"
            f"<td>{html.escape(interp)}</td><td class='num'>{cnt:,}</td>"
            f"<td class='num'>{pct_lib:.2f}%</td><td class='num'>{pct_all:.2f}%</td></tr>"
        )

    sort_js = """<script>
document.querySelectorAll('table').forEach(tbl=>{
  const ths=tbl.querySelectorAll('thead th');
  ths.forEach((th,idx)=>{
    th.dataset.orig=th.textContent.trim();
    th.onclick=()=>{
      const tb=tbl.querySelector('tbody'); if(!tb)return;
      const asc=th.dataset.dir!=='asc';
      ths.forEach(h=>{h.dataset.dir='';if(h.dataset.orig)h.textContent=h.dataset.orig;});
      th.dataset.dir=asc?'asc':'desc';
      th.textContent=th.dataset.orig+(asc?' ▲':' ▼');
      const pv=td=>{const s=(td?td.textContent:'').trim(),c=s.replace(/[,%\\sx×]/gi,'');return /^-?\\d+(\\.\\d+)?$/.test(c)&&c!==''?{n:parseFloat(c),ok:true,s}:{n:0,ok:false,s};};
      const rows=Array.from(tb.querySelectorAll('tr'));
      rows.sort((a,b)=>{const va=pv(a.children[idx]),vb=pv(b.children[idx]);return(va.ok&&vb.ok)?(asc?va.n-vb.n:vb.n-va.n):(asc?va.s.localeCompare(vb.s,'zh-Hant'):vb.s.localeCompare(va.s,'zh-Hant'));});
      rows.forEach(r=>tb.appendChild(r));
    };
  });
});
</script>"""

    sankey_html = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<title>{html.escape(sample_id)} — 全序列最終判定 Sankey 圖與明細表</title>
<style>{base_css}</style>
</head>
<body>
<div class="container">
  <div class="header-card">
    <div>
      <h1>ClinMeta <code>{html.escape(sample_id)}</code> — 全序列最終判定 Sankey 圖與明細表</h1>
      <p class="sub">完整追蹤全數 <strong>{total_raw:,}</strong> 條原始定序 reads（{lib_desc}）從 QC、去宿主到最終病原體／背景／人為產物判定的互斥歸屬。</p>
    </div>
    {logo_html}
  </div>
  <div class="panel">
    <div style="display:flex;flex-wrap:wrap;gap:12px;align-items:center;margin-bottom:12px;padding:10px 12px;background:#f8fafc;border:1px solid #e2e8f0;border-radius:8px;font-size:13px;">
      <label>顯示範圍：
        <select id="sk-scope" onchange="renderInteractiveSankey()" style="padding:4px 8px;border:1px solid #cbd5e1;border-radius:6px;">
          <option value="all">全部序列（含 QC 與宿主）</option>
          <option value="nonhost">僅看非人類（Non-Host 流向）</option>
        </select>
      </label>
      <label>比例縮放：
        <select id="sk-scale" onchange="renderInteractiveSankey()" style="padding:4px 8px;border:1px solid #cbd5e1;border-radius:6px;">
          <option value="pow">冪次縮放（兼顧微量病原）</option>
          <option value="log">對數縮放（Log10）</option>
          <option value="linear">線性真實比例</option>
        </select>
      </label>
      <label>標籤顯示：
        <select id="sk-lbl" onchange="renderInteractiveSankey()" style="padding:4px 8px;border:1px solid #cbd5e1;border-radius:6px;">
          <option value="count">Reads 數</option>
          <option value="pct">佔總 Raw 百分比</option>
          <option value="both">Reads + 百分比</option>
        </select>
      </label>
      <button onclick="downloadSankeySvg()" style="margin-left:auto;background:#2563eb;color:#fff;border:none;border-radius:6px;padding:6px 12px;font-size:12.5px;font-weight:600;cursor:pointer;">匯出 SVG</button>
    </div>
    <div id="sankey-box">{sankey_svg}</div>
  </div>
  <script>
    const SK_RAW_NODES = {json.dumps(nodes, ensure_ascii=False)};
    const SK_RAW_LINKS = {json.dumps(links, ensure_ascii=False)};
    const SK_TOTAL = {max(1, total_raw)};
    const SK_SID = {json.dumps(sample_id, ensure_ascii=False)};
    function renderInteractiveSankey() {{
      const scope = document.getElementById('sk-scope').value;
      const scaleMode = document.getElementById('sk-scale').value;
      const lblMode = document.getElementById('sk-lbl').value;
      const nodes = JSON.parse(JSON.stringify(SK_RAW_NODES));
      let links = JSON.parse(JSON.stringify(SK_RAW_LINKS));
      if (scope === 'nonhost') {{
        links = links.filter(l => l.source === 'nonhost');
        const keep = new Set(['nonhost', ...links.map(l => l.target)]);
        Object.keys(nodes).forEach(k => {{ if (!keep.has(k)) delete nodes[k]; else if (nodes[k].col >= 2) nodes[k].col -= 2; }});
      }}
      const W = 1160, H = 520, nodeW = 18, padTop = 46, padBot = 22, availH = H - padTop - padBot;
      const cols = [...new Set(Object.values(nodes).map(n => n.col))].sort((a,b) => a-b);
      if (!cols.length) return;
      const colX = {{}};
      cols.forEach((c, i) => {{ colX[c] = 24 + i * ((W - 360) / Math.max(1, cols.length - 1)); }});
      Object.entries(nodes).forEach(([id, n]) => {{
        const inV = links.filter(l => l.target === id).reduce((s, l) => s + l.val, 0);
        const outV = links.filter(l => l.source === id).reduce((s, l) => s + l.val, 0);
        n.val = Math.max(inV, outV, 0);
      }});
      const sv = v => scaleMode === 'linear' ? Math.max(v, 1) : (scaleMode === 'log' ? Math.log10(Math.max(v, 1) + 1) : Math.pow(Math.max(v, 1), 0.45));
      cols.forEach(c => {{
        const cnodes = Object.keys(nodes).filter(id => nodes[id].col === c);
        const gap = 12, totS = cnodes.reduce((s, id) => s + sv(nodes[id].val), 0) || 1;
        const usableH = Math.max(40, availH - gap * Math.max(0, cnodes.length - 1));
        let y = padTop;
        cnodes.forEach(id => {{
          const nh = Math.max(6, usableH * (sv(nodes[id].val) / totS));
          Object.assign(nodes[id], {{ x: colX[c], y, h: nh, outY: y, inY: y }});
          y += nh + gap;
        }});
      }});
      const outT = {{}}, inT = {{}};
      Object.keys(nodes).forEach(id => {{
        outT[id] = links.filter(l => l.source === id).reduce((s, l) => s + l.val, 0) || 1;
        inT[id] = links.filter(l => l.target === id).reduce((s, l) => s + l.val, 0) || 1;
      }});
      const esc = s => String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
      const paths = links.map(l => {{
        const s = nodes[l.source], t = nodes[l.target];
        const sh = s.h * (l.val / outT[l.source]), th = t.h * (l.val / inT[l.target]);
        const y0 = s.outY + sh/2, y1 = t.inY + th/2;
        s.outY += sh; t.inY += th;
        const x0 = s.x + nodeW, x1 = t.x, xm = (x0 + x1) / 2, sw = Math.max(2, (sh + th) / 2);
        const tip = `${{s.label}} -> ${{t.label}}: ${{l.val.toLocaleString()}} reads (${{(l.val/SK_TOTAL*100).toFixed(3)}}%)`;
        return `<path d="M ${{x0.toFixed(1)}},${{y0.toFixed(1)}} C ${{xm.toFixed(1)}},${{y0.toFixed(1)}} ${{xm.toFixed(1)}},${{y1.toFixed(1)}} ${{x1.toFixed(1)}},${{y1.toFixed(1)}}" fill="none" stroke="${{l.color||'#94a3b8'}}" stroke-opacity="0.38" stroke-width="${{sw.toFixed(1)}}" class="sankey-link"><title>${{esc(tip)}}</title></path>`;
      }}).join('');
      const rects = Object.values(nodes).map(n => {{
        const pct = (n.val / SK_TOTAL * 100).toFixed(2) + '%';
        const valTxt = lblMode === 'pct' ? pct : (lblMode === 'both' ? `${{n.val.toLocaleString()}} (${{pct}})` : n.val.toLocaleString());
        const lbl = `${{n.label}} (${{valTxt}})`;
        return `<rect x="${{n.x.toFixed(1)}}" y="${{n.y.toFixed(1)}}" width="${{nodeW}}" height="${{n.h.toFixed(1)}}" rx="3" fill="${{n.color||'#334155'}}"><title>${{esc(lbl)}}</title></rect><text x="${{(n.x+nodeW+6).toFixed(1)}}" y="${{(n.y+n.h/2+4).toFixed(1)}}" font-size="11.5" font-weight="600" fill="#0f172a">${{esc(lbl)}}</text>`;
      }}).join('');
      const title = scope === 'nonhost' ? `${{SK_SID}} 非人類 (Non-Host) 序列流向 Sankey 圖` : `${{SK_SID}} 全序列流向與最終判定 Sankey 圖（總計 ${{SK_TOTAL.toLocaleString()}} Raw Reads）`;
      document.getElementById('sankey-box').innerHTML = `<svg viewBox="0 0 ${{W}} ${{H}}" style="width:100%;height:auto;background:#fff;border:1px solid #e2e8f0;border-radius:10px;"><text x="20" y="26" font-size="14" font-weight="700" fill="#0f172a">${{esc(title)}}</text>${{paths}}${{rects}}</svg>`;
    }}
    function downloadSankeySvg() {{
      const svg = document.querySelector('#sankey-box svg');
      if (!svg) return;
      const blob = new Blob([svg.outerHTML], {{type:'image/svg+xml;charset=utf-8'}});
      const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = SK_SID + '_sankey.svg'; a.click();
    }}
    renderInteractiveSankey();
  </script>
  <div class="panel">
    <h2>全序列互斥最終判定明細表（合計 {total_raw:,} reads = 100.00%）</h2>
    <table>
      <thead><tr><th>Library</th><th>最終分類項目</th><th>判定方法</th><th>生物／臨床意義</th><th class="num">Reads</th><th class="num">佔該 Library</th><th class="num">佔總 Raw</th></tr></thead>
      <tbody>{''.join(d_html_rows)}</tbody>
    </table>
  </div>
</div>
{sort_js}
</body>
</html>"""
    with open(os.path.join(outdir, "read_disposition_sankey.html"), "w") as f:
        f.write(sankey_html)

    # --- 2. Tab 2: Pathogen Competitive Mapping Validation Report ---
    ranked = sorted(species_rows, key=lambda r: (r["mapped_reads_mapq20"] >= 10, r["breadth_pct"], r["mapped_reads_mapq20"]), reverse=True)
    confirmed = [r for r in ranked if r["mapped_reads_mapq20"] >= 10]
    unique_confirmed = []
    for r in confirmed:
        if r["species"] not in [x["species"] for x in unique_confirmed]:
            unique_confirmed.append(r)

    if unique_confirmed:
        p1 = unique_confirmed[0]
        card1 = (
            f"<div class='card green'><div class='card-title'>主要檢出病原體</div>"
            f"<div class='card-val'>{html.escape(p1['species'])}</div>"
            f"<div class='card-note'>{html.escape(p1['library'])} {p1['mapped_reads']:,} reads（MAPQ≥20: {p1['mapped_reads_mapq20']:,}），"
            f"覆蓋廣度 <strong>{p1['breadth_pct']:.2f}%</strong>，平均深度 <strong>{p1['mean_depth']:.2f}×</strong>。</div></div>"
        )
    else:
        card1 = (
            "<div class='card slate'><div class='card-title'>主要檢出病原體</div>"
            "<div class='card-val'>未檢出顯著病原體</div>"
            "<div class='card-note'>所有候選參考基因組之 MAPQ≥20 reads 皆低於門檻（&lt;10 reads）。</div></div>"
        )

    if len(unique_confirmed) > 1:
        p2 = unique_confirmed[1]
        card2 = (
            f"<div class='card green'><div class='card-title'>次要檢出／共感染訊號</div>"
            f"<div class='card-val'>{html.escape(p2['species'])}</div>"
            f"<div class='card-note'>{html.escape(p2['library'])} {p2['mapped_reads']:,} reads（MAPQ≥20: {p2['mapped_reads_mapq20']:,}），"
            f"覆蓋廣度 <strong>{p2['breadth_pct']:.2f}%</strong>，平均深度 <strong>{p2['mean_depth']:.2f}×</strong>。</div></div>"
        )
    else:
        card2 = (
            "<div class='card amber'><div class='card-title'>次要／背景交叉比對</div>"
            f"<div class='card-val'>共評估 {len(ranked)} 筆物種</div>"
            "<div class='card-note'>經競爭性比對索引過濾後，無其他顯著次要共感染病原體。</div></div>"
        )

    host_pct = (tot_host / max(1, total_raw) * 100.0)
    card3 = (
        f"<div class='card'><div class='card-title'>去宿主與非宿主保留</div>"
        f"<div class='card-val'>{tot_nonhost:,} non-host</div>"
        f"<div class='card-note'>原始 {total_raw:,} reads；Hostile 去除人類宿主 {tot_host:,}（{host_pct:.1f}%），QC 淘汰 {tot_qc:,}。</div></div>"
    )
    card4 = (
        f"<div class='card slate'><div class='card-title'>競爭性回貼分流</div>"
        f"<div class='card-val'>{tot_pathogen:,} pathogen reads</div>"
        f"<div class='card-note'>殘餘人類回貼 {tot_res_human:,}；未比對至病原體面板 {tot_unmapped:,}。</div></div>"
    )

    cov_bars = []
    for r in ranked:
        q20 = int(r["mapped_reads_mapq20"])
        cov = float(r["breadth_pct"])
        bar_cls = "bar-fill" if q20 >= 10 else "bar-fill warn"
        w_pct = min(100.0, max(2.0 if cov > 0 else 0.0, cov))
        cov_bars.append(
            f"<div class='bar-row'>"
            f"<span class='bar-label' title='{html.escape(r['library'])} | {html.escape(r['species'])}'>[{html.escape(r['library'])}] {html.escape(r['species'])}</span>"
            f"<div class='bar-track'><div class='{bar_cls}' style='width:{w_pct:.1f}%'></div></div>"
            f"<span class='bar-val'>{cov:.2f}% cov · {q20:,} Q20 · {float(r['mean_depth']):.2f}×</span>"
            f"</div>"
        )

    sp_html_rows = []
    for r in ranked:
        q20 = int(r["mapped_reads_mapq20"])
        cov = float(r["breadth_pct"])
        if q20 >= 10 and cov >= 1.0:
            badge = "<span class='badge badge-pos'>Confirmed</span>"
        elif q20 >= 10:
            badge = "<span class='badge badge-pos'>Detected</span>"
        elif q20 > 0:
            badge = "<span class='badge badge-warn'>Low / Trace</span>"
        else:
            badge = "<span class='badge badge-neg'>Not Detected</span>"
        sp_html_rows.append(
            f"<tr><td>{html.escape(r['library'])}</td>"
            f"<td><strong>{html.escape(r['species'])}</strong>{badge}</td>"
            f"<td class='num'>{int(r['ref_len']):,}</td>"
            f"<td class='num'>{int(r['mapped_reads']):,}</td>"
            f"<td class='num'>{q20:,}</td>"
            f"<td class='num'>{cov:.2f}%</td>"
            f"<td class='num'>{float(r['mean_depth']):.2f}x</td></tr>"
        )

    val_html = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<title>{html.escape(sample_id)} — 病原體競爭性回貼驗證報告</title>
<style>{base_css}</style>
</head>
<body>
<div class="container">
  <div class="header-card">
    <div>
      <h1>ClinMeta <code>{html.escape(sample_id)}</code> — 病原體競爭性回貼驗證報告</h1>
      <p class="sub">競爭性比對索引：<code>human-masked.fa.gz</code> + Offline Pathogen Panel 候選參考基因組（<code>minimap2 -ax map-ont --secondary=no</code>）</p>
    </div>
    {logo_html}
  </div>
  <div class="cards">{card1}{card2}{card3}{card4}</div>
  <div class="panel">
    <h2>物種基因組覆蓋廣度與高可信讀段比較（MAPQ≥20）</h2>
    {''.join(cov_bars) if cov_bars else '<p class="sub">無候選參考物種進行競爭性回貼。</p>'}
  </div>
  <div class="panel">
    <h2>離線競爭性回貼驗證明細表</h2>
    <table>
      <thead><tr><th>Library</th><th>物種與判定</th><th class="num">參考長度</th><th class="num">Mapped Reads</th><th class="num">MAPQ≥20</th><th class="num">覆蓋廣度</th><th class="num">平均深度</th></tr></thead>
      <tbody>{''.join(sp_html_rows) or '<tr><td colspan="7">無顯著匹配病原體</td></tr>'}</tbody>
    </table>
  </div>
</div>
{sort_js}
</body>
</html>"""
    with open(os.path.join(outdir, "pathogen_validation_report.html"), "w") as f:
        f.write(val_html)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample-id", required=True)
    p.add_argument("--dna-path", default="")
    p.add_argument("--cdna-path", default="")
    p.add_argument("--specimen-type", default="unknown")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--db-dir", default="databases")
    p.add_argument("--ref-dir", default="references")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    dna_fqs = find_fastqs(args.dna_path)
    cdna_fqs = find_fastqs(args.cdna_path)
    if not dna_fqs and not cdna_fqs:
        raise SystemExit("Error: no DNA or cDNA FASTQ files found.")

    os.makedirs(args.outdir, exist_ok=True)
    dehost_dir = os.path.join(args.outdir, "dehost")
    k2_dir = os.path.join(args.outdir, "kraken2")
    tmp_dir = os.path.join(args.outdir, "tmp")
    os.makedirs(dehost_dir, exist_ok=True)
    os.makedirs(k2_dir, exist_ok=True)

    plan = compute_execution_plan(queued_count=1, requested_threads=args.threads)
    threads = args.threads if args.threads > 0 else plan["threads_per_job"]

    # 1. QC & Dehost (run DNA and cDNA in parallel when both exist)
    update_progress(args.outdir, "1/4 QC & Dehost", 15, f"Running chopper QC and hostile dehosting ({threads} threads)...")
    dna_clean = os.path.join(dehost_dir, f"{args.sample_id}.DNA.nonhost.fastq.gz")
    cdna_clean = os.path.join(dehost_dir, f"{args.sample_id}.cDNA.nonhost.fastq.gz")
    dehost_stats = {}
    libs_clean = []

    tasks = []
    if dna_fqs:
        tasks.append(("DNA", dna_fqs, "map-ont", os.path.join(args.ref_dir, "human-masked.mmi"), dna_clean))
    if cdna_fqs:
        tasks.append(("cDNA", cdna_fqs, "splice", os.path.join(args.ref_dir, "human-masked.splice.mmi"), cdna_clean))
    per_lib_threads = max(4, threads // max(1, len(tasks)))
    with ThreadPoolExecutor(max_workers=max(1, len(tasks))) as pool:
        futs = {
            lbl: pool.submit(run_qc_and_dehost, fqs, lbl, preset, idx, out_fq, per_lib_threads, tmp_dir)
            for lbl, fqs, preset, idx, out_fq in tasks
        }
        for lbl, _, _, _, out_fq in tasks:
            dehost_stats[lbl] = futs[lbl].result()
            libs_clean.append((lbl, out_fq))

    # 2. Kraken2
    update_progress(args.outdir, "2/4 Kraken2", 45, "Running Kraken2 PlusPF classification...")
    k2_db = os.path.join(args.db_dir, "kraken2", "k2_pluspf_20261022")
    k2_reports = []
    for label, fq in libs_clean:
        rep = os.path.join(k2_dir, f"{args.sample_id}.{label}.report.txt")
        out_k = os.path.join(k2_dir, f"{args.sample_id}.{label}.reads.kraken")
        if count_fastq_reads(fq) == 0:
            open(rep, "w").close()
            with gzip.open(out_k + ".gz", "wb"):
                pass
        else:
            k2_cmd = [
                "kraken2", "--db", k2_db, "--threads", str(threads),
                "--report-minimizer-data", "--gzip-compressed",
                "--report", rep, "--output", out_k,
            ]
            if plan["use_mmap_kraken"]:
                k2_cmd.append("--memory-mapping")
            k2_cmd.append(fq)
            subprocess.run(k2_cmd, check=True)
            if os.path.exists(out_k):
                subprocess.run(["pigz", "-f", "-p", "4", out_k], check=True)
            else:
                with gzip.open(out_k + ".gz", "wb"):
                    pass
        k2_reports.append(rep)

    # Generate interactive Plotly/Sankey Kraken2 HTML report
    k2_html_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kraken2_report_html.py")
    if os.path.exists(k2_html_script) and k2_reports:
        subprocess.run([sys.executable, k2_html_script, k2_dir, os.path.join(k2_dir, "report.html")], check=False)

    # 3. Offline competitive validation
    update_progress(args.outdir, "3/4 Competitive Validation", 75, "Running offline competitive reference mapping...")
    panel_dir = os.path.join(args.ref_dir, "pathogen_panel")
    selected_refs = select_panel_refs(k2_reports, panel_dir)
    species_rows, map_overview = run_competitive_validation(
        libs_clean, selected_refs, os.path.join(args.ref_dir, "human-masked.fa.gz"),
        args.outdir, threads, use_shm_tmp=plan["use_shm_tmp"],
    )

    # Determine primary finding summary
    top_findings = sorted(
        [r for r in species_rows if r["mapped_reads_mapq20"] >= 10],
        key=lambda r: (r["breadth_pct"], r["mapped_reads_mapq20"]),
        reverse=True,
    )
    seen_sp = []
    for r in top_findings:
        if r["species"] not in [x[0] for x in seen_sp]:
            seen_sp.append((r["species"], f"{r['species']} ({r['library']} {r['breadth_pct']:.1f}% cov, {r['mapped_reads_mapq20']} Q20)"))
        if len(seen_sp) >= 2:
            break
    primary_finding = " + ".join(x[1] for x in seen_sp) if seen_sp else "No major pathogen above threshold"

    # 4. Reports & Finalize
    update_progress(args.outdir, "4/4 Reports", 92, "Generating reports...")
    write_simple_reports(args.sample_id, dehost_stats, species_rows, map_overview, args.outdir)

    summary = {
        "sample_id": args.sample_id,
        "specimen_type": args.specimen_type,
        "primary_finding": primary_finding,
        "dehost_stats": dehost_stats,
    }
    with open(os.path.join(args.outdir, "sample_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    update_progress(args.outdir, "Completed", 100, "Analysis complete.", summary)


if __name__ == "__main__":
    main()
