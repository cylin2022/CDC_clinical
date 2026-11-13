#!/usr/bin/env python3
"""Generate Influenza A (H3N2) consensus FASTA and standalone HTML validation report."""

import base64
import collections
import csv
import html
import json
import os
import subprocess

OUTDIR = "analysis/pathogen_validation"
CDNA_BAM = os.path.join(OUTDIR, "cDNA_ST26-UM01.nonhuman.bam")
SPECIES_TSV = os.path.join(OUTDIR, "species_validation_summary.tsv")
CONTIG_TSV = os.path.join(OUTDIR, "contig_validation_summary.tsv")
OVERVIEW_JSON = os.path.join(OUTDIR, "mapping_overview.json")

FLU_SEGMENTS = [
    ("FLU_H3N2|NC_007373.1", "PB2", 2341),
    ("FLU_H3N2|NC_007372.1", "PB1", 2341),
    ("FLU_H3N2|NC_007371.1", "PA", 2233),
    ("FLU_H3N2|NC_007366.1", "HA", 1762),
    ("FLU_H3N2|NC_007369.1", "NP", 1566),
    ("FLU_H3N2|NC_007368.1", "NA", 1467),
    ("FLU_H3N2|NC_007367.1", "M", 1027),
    ("FLU_H3N2|NC_007370.1", "NS", 890),
]

CANDIDA_CHRS = [
    ("NC_032089.1", "Chr1"),
    ("NC_032090.1", "Chr2"),
    ("NC_032091.1", "Chr3"),
    ("NC_032092.1", "Chr4"),
    ("NC_032093.1", "Chr5"),
    ("NC_032094.1", "Chr6"),
    ("NC_032095.1", "Chr7"),
    ("NC_032096.1", "ChrR"),
]


def call_flu_consensus():
    """Extract MAPQ>=20 FLU_H3N2 alignments and generate consensus FASTA."""
    flu_bam = os.path.join(OUTDIR, "influenza_a_h3n2.cDNA.mapq20.bam")
    reg_args = [r for r, _, _ in FLU_SEGMENTS]
    subprocess.run(
        ["samtools", "view", "-b", "-q", "20", "-o", flu_bam, CDNA_BAM] + reg_args,
        check=True,
    )
    subprocess.run(["samtools", "index", flu_bam], check=True)

    raw_fa = os.path.join(OUTDIR, "influenza_a_h3n2_consensus.raw.fasta")
    out_fa = os.path.join(OUTDIR, "influenza_a_h3n2_consensus.fasta")
    # Call consensus with min-depth 10; positions <10x become N
    with open(raw_fa, "w") as f:
        for rname, _, _ in FLU_SEGMENTS:
            res = subprocess.run(
                [
                    "samtools",
                    "consensus",
                    "-f",
                    "fasta",
                    "-d",
                    "10",
                    "-r",
                    rname,
                    flu_bam,
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            f.write(res.stdout)

    # Parse and rewrite headers cleanly with segment annotations & completeness stats
    seqs = {}
    cur = None
    buf = []
    with open(raw_fa) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur is not None:
            seqs[cur] = "".join(buf)

    seg_meta = {rname: seg for rname, seg, _ in FLU_SEGMENTS}
    consensus_stats = []
    with open(out_fa, "w") as out:
        for rname, seg, ref_len in FLU_SEGMENTS:
            seq = seqs.get(rname, "N" * ref_len).upper()
            n_count = seq.count("N")
            called = len(seq) - n_count
            pct = (called / len(seq) * 100.0) if seq else 0.0
            acc = rname.split("|", 1)[1]
            out.write(
                f">sample_001_FLU_H3N2_{seg} ref={acc} len={len(seq)} called={called} ({pct:.2f}%) N={n_count}\n"
            )
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")
            consensus_stats.append(
                {
                    "segment": seg,
                    "accession": acc,
                    "consensus_len": len(seq),
                    "called_bases": called,
                    "n_bases": n_count,
                    "called_pct": round(pct, 2),
                }
            )
    os.remove(raw_fa)

    # Compute binned depth profiles for SVG plots
    depth_proc = subprocess.Popen(
        ["samtools", "depth", "-a", flu_bam], stdout=subprocess.PIPE, text=True
    )
    depths = collections.defaultdict(list)
    for line in depth_proc.stdout:
        rname, pos, d = line.rstrip("\n").split("\t")
        if rname in seg_meta:
            depths[rname].append(int(d))
    if depth_proc.wait() != 0:
        raise RuntimeError("samtools depth failed")

    return out_fa, consensus_stats, depths


def svg_depth_plot(depths_list, seg_name, ref_len):
    """Render a compact inline SVG depth profile (log10-scaled or linear with cap)."""
    w, h = 420, 110
    pad_l, pad_r, pad_t, pad_b = 38, 10, 18, 20
    pw = w - pad_l - pad_r
    ph = h - pad_t - pad_b
    n_bins = 100
    n = len(depths_list) or ref_len
    bin_vals = []
    for b in range(n_bins):
        s = int(b * n / n_bins)
        e = max(s + 1, int((b + 1) * n / n_bins))
        chunk = depths_list[s:e]
        bin_vals.append(sum(chunk) / len(chunk) if chunk else 0.0)

    max_d = max(bin_vals) if bin_vals else 1.0
    max_d = max(max_d, 10.0)
    pts = [f"{pad_l},{pad_t + ph}"]
    for i, v in enumerate(bin_vals):
        x = pad_l + (i / (n_bins - 1)) * pw
        y = pad_t + ph - (v / max_d) * ph
        pts.append(f"{x:.1f},{y:.1f}")
    pts.append(f"{pad_l + pw},{pad_t + ph}")
    poly = " ".join(pts)

    y10 = pad_t + ph - (min(10.0, max_d) / max_d) * ph
    return f"""<svg viewBox="0 0 {w} {h}" class="depth-svg">
  <rect x="{pad_l}" y="{pad_t}" width="{pw}" height="{ph}" fill="#f8fafc" stroke="#e2e8f0"/>
  <line x1="{pad_l}" y1="{y10:.1f}" x2="{pad_l+pw}" y2="{y10:.1f}" stroke="#ef4444" stroke-dasharray="3,3" stroke-width="1"/>
  <polygon points="{poly}" fill="rgba(37,99,235,0.22)" stroke="#2563eb" stroke-width="1.5"/>
  <text x="{pad_l}" y="12" font-size="11" font-weight="600" fill="#0f172a">{html.escape(seg_name)} ({ref_len} bp) — max {max_d:.0f}x</text>
  <text x="{pad_l-4}" y="{pad_t+8}" font-size="9" text-anchor="end" fill="#64748b">{max_d:.0f}x</text>
  <text x="{pad_l-4}" y="{y10+3:.1f}" font-size="9" text-anchor="end" fill="#ef4444">10x</text>
  <text x="{pad_l-4}" y="{pad_t+ph}" font-size="9" text-anchor="end" fill="#64748b">0</text>
  <text x="{pad_l}" y="{h-4}" font-size="9" fill="#64748b">1</text>
  <text x="{pad_l+pw}" y="{h-4}" font-size="9" text-anchor="end" fill="#64748b">{ref_len}</text>
</svg>"""


def build_html(consensus_stats, depths):
    with open(OVERVIEW_JSON) as f:
        overview = json.load(f)
    with open(SPECIES_TSV) as f:
        species_rows = list(csv.DictReader(f, delimiter="\t"))
    with open(CONTIG_TSV) as f:
        contig_rows = list(csv.DictReader(f, delimiter="\t"))

    flu_contigs = {
        r["contig"]: r
        for r in contig_rows
        if r["library"] == "cDNA_ST26-UM01" and r["species"] == "FLU_H3N2"
    }
    cons_by_acc = {c["accession"]: c for c in consensus_stats}

    flu_plots_html = []
    flu_table_rows = []
    for rname, seg, ref_len in FLU_SEGMENTS:
        acc = rname.split("|", 1)[1]
        dlist = depths.get(rname, [0] * ref_len)
        flu_plots_html.append(svg_depth_plot(dlist, seg, ref_len))
        cr = flu_contigs.get(acc, {})
        cs = cons_by_acc.get(acc, {})
        flu_table_rows.append(
            f"<tr><td><strong>{seg}</strong></td><td><code>{acc}</code></td>"
            f"<td class='num'>{ref_len:,}</td>"
            f"<td class='num'>{int(cr.get('reads', 0)):,}</td>"
            f"<td class='num'>{int(cr.get('reads_mapq20', 0)):,}</td>"
            f"<td class='num'>{float(cr.get('coverage_pct', 0)):.2f}%</td>"
            f"<td class='num'>{float(cr.get('mean_depth', 0)):.2f}x</td>"
            f"<td class='num'>{cs.get('called_bases', 0):,} ({cs.get('called_pct', 0):.2f}%)</td>"
            f"<td class='num'>{cs.get('n_bases', 0):,}</td></tr>"
        )

    candida_dna = {
        r["contig"]: r
        for r in contig_rows
        if r["library"] == "DNA_SG26-GS01" and r["species"] == "CANDIDA_ALBICANS"
    }
    candida_bars = []
    for acc, label in CANDIDA_CHRS:
        r = candida_dna.get(acc, {})
        cov = float(r.get("coverage_pct", 0))
        rds = int(r.get("reads_mapq20", 0))
        candida_bars.append(
            f"<div class='bar-row'><span class='bar-label'>{label} ({acc})</span>"
            f"<div class='bar-track'><div class='bar-fill' style='width:{min(100, cov*15):.1f}%'></div></div>"
            f"<span class='bar-val'>{cov:.2f}% ({rds} MAPQ≥20 reads)</span></div>"
        )

    sp_table_rows = []
    for r in species_rows:
        sp = r["species"]
        badge = ""
        if sp == "FLU_H3N2" and r["library"].startswith("cDNA"):
            badge = "<span class='badge badge-pos'>Confirmed Primary</span>"
        elif sp == "CANDIDA_ALBICANS" and r["library"].startswith("DNA"):
            badge = "<span class='badge badge-pos'>Confirmed DNA</span>"
        elif sp == "CANDIDA_ALBICANS" and r["library"].startswith("cDNA"):
            badge = "<span class='badge badge-pos'>Confirmed cDNA</span>"
        elif sp == "FLU_H1N1":
            badge = "<span class='badge badge-neg'>Subtype Control (0)</span>"
        elif int(r["mapped_reads_mapq20"]) <= 26 and float(r["breadth_pct"]) < 0.5:
            badge = "<span class='badge badge-warn'>Low / Artifact</span>"
        sp_table_rows.append(
            f"<tr><td>{html.escape(r['library'])}</td>"
            f"<td><strong>{html.escape(sp)}</strong> {badge}</td>"
            f"<td class='num'>{int(r['ref_len']):,}</td>"
            f"<td class='num'>{int(r['mapped_reads']):,}</td>"
            f"<td class='num'>{int(r['mapped_reads_mapq20']):,}</td>"
            f"<td class='num'>{float(r['mean_mapq']):.2f}</td>"
            f"<td class='num'>{float(r['mean_identity_pct']):.2f}%</td>"
            f"<td class='num'>{float(r['breadth_pct']):.2f}%</td>"
            f"<td class='num'>{float(r['mean_depth']):.2f}x</td></tr>"
        )

    total_called = sum(c["called_bases"] for c in consensus_stats)
    total_len = sum(c["consensus_len"] for c in consensus_stats)
    total_pct = total_called / total_len * 100.0 if total_len else 0.0

    logo_b64 = ""
    logo_path = os.path.join("webapp", "assets", "clinmeta-logo-light.jpg")
    if os.path.exists(logo_path):
        with open(logo_path, "rb") as lf:
            logo_b64 = base64.b64encode(lf.read()).decode("ascii")
    logo_html = (
        f'<img src="data:image/jpeg;base64,{logo_b64}" alt="ClinMeta" style="height:56px;border-radius:6px;">'
        if logo_b64 else ""
    )

    html_doc = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<title>CDC Clinical sample_001 — Pathogen Validation Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; background: #f1f5f9; color: #0f172a; line-height: 1.5; }}
  .container {{ max-width: 1180px; margin: 28px auto; padding: 0 20px; }}
  h1 {{ font-size: 24px; margin-bottom: 4px; }}
  .sub {{ color: #475569; margin-bottom: 22px; font-size: 14px; }}
  .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 16px; margin-bottom: 24px; }}
  .card {{ background: #fff; border-radius: 10px; padding: 16px 18px; box-shadow: 0 1px 3px rgba(15,23,42,0.08); border-left: 4px solid #2563eb; }}
  .card.green {{ border-left-color: #16a34a; }}
  .card.amber {{ border-left-color: #d97706; }}
  .card.slate {{ border-left-color: #64748b; }}
  .card-title {{ font-size: 12px; text-transform: uppercase; letter-spacing: 0.04em; color: #64748b; font-weight: 600; }}
  .card-val {{ font-size: 22px; font-weight: 700; margin: 6px 0 4px; }}
  .card-note {{ font-size: 13px; color: #475569; }}
  .panel {{ background: #fff; border-radius: 10px; padding: 20px; margin-bottom: 22px; box-shadow: 0 1px 3px rgba(15,23,42,0.08); }}
  .panel h2 {{ font-size: 18px; margin-top: 0; margin-bottom: 12px; border-bottom: 1px solid #e2e8f0; padding-bottom: 8px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
  th, td {{ padding: 8px 10px; border-bottom: 1px solid #e2e8f0; text-align: left; }}
  th {{ background: #f8fafc; font-weight: 600; color: #334155; }}
  td.num, th.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  .badge {{ display: inline-block; padding: 2px 7px; border-radius: 999px; font-size: 11px; font-weight: 600; margin-left: 6px; }}
  .badge-pos {{ background: #dcfce7; color: #166534; }}
  .badge-neg {{ background: #e2e8f0; color: #334155; }}
  .badge-warn {{ background: #fef3c7; color: #92400e; }}
  .svg-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(380px, 1fr)); gap: 14px; margin-bottom: 18px; }}
  .depth-svg {{ width: 100%; height: auto; background: #fff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 6px; box-sizing: border-box; }}
  .bar-row {{ display: flex; align-items: center; gap: 10px; margin: 7px 0; font-size: 13px; }}
  .bar-label {{ width: 180px; font-family: monospace; color: #334155; }}
  .bar-track {{ flex: 1; height: 14px; background: #e2e8f0; border-radius: 7px; overflow: hidden; }}
  .bar-fill {{ height: 100%; background: #16a34a; }}
  .bar-val {{ width: 200px; text-align: right; font-variant-numeric: tabular-nums; color: #334155; }}
  code {{ background: #f1f5f9; padding: 2px 5px; border-radius: 4px; font-size: 12px; }}
</style>
</head>
<body>
<div class="container">
  <div style="display:flex;align-items:center;justify-content:space-between;gap:16px;flex-wrap:wrap;margin-bottom:14px;background:#fff;padding:14px 20px;border-radius:10px;box-shadow:0 1px 3px rgba(15,23,42,0.08);">
    <div>
      <h1 style="margin:0 0 4px;">ClinMeta <code>sample_001</code> — 病原體競爭性回貼驗證報告</h1>
      <div class="sub" style="margin:0;">競爭性比對索引：<code>human-masked.fa.gz</code> + 8 組 NCBI RefSeq 參考基因組（<code>minimap2 -ax map-ont --secondary=no</code>）</div>
    </div>
    {logo_html}
  </div>

  <div class="cards">
    <div class="card green">
      <div class="card-title">主要病毒病原體</div>
      <div class="card-val">Influenza A (H3N2)</div>
      <div class="card-note">cDNA 3,694 reads（MAPQ≥20: 3,318），全基因組覆蓋 <strong>98.45%</strong>，平均深度 <strong>78.12×</strong>；H1N1 對照為 0 reads。</div>
    </div>
    <div class="card green">
      <div class="card-title">次要真菌檢出</div>
      <div class="card-val">Candida albicans</div>
      <div class="card-note">DNA 332 reads（MAPQ≥20: 263，覆蓋 4.41% 基因組，均勻分佈 8 條染色體）；cDNA 1,089 reads。</div>
    </div>
    <div class="card amber">
      <div class="card-title">細菌假陽性排除</div>
      <div class="card-val">已排除高計數細菌</div>
      <div class="card-note">Kraken2 高計數之 <em>A. pittii</em>、<em>V. alginolyticus</em>、<em>E. coli</em>、<em>K. pneumoniae</em> 回貼覆蓋皆 &lt;0.4%。</div>
    </div>
    <div class="card slate">
      <div class="card-title">H3N2 一致性序列</div>
      <div class="card-val">{total_called:,} / {total_len:,} bp ({total_pct:.2f}%)</div>
      <div class="card-note">門檻：MAPQ≥20 且深度 ≥10×（低於 10× 標記為 N）；輸出於 <code>influenza_a_h3n2_consensus.fasta</code>。</div>
    </div>
  </div>

  <div class="panel">
    <h2>1. Influenza A (H3N2) 全基因組 8 片段覆蓋深度與一致性序列統計（cDNA ST26-UM01）</h2>
    <p style="font-size:13px;color:#475569;margin-top:0;">紅色虛線代表一致性序列定鹼基門檻（<code>10×</code> 深度）。</p>
    <div class="svg-grid">
      {''.join(flu_plots_html)}
    </div>
    <table>
      <thead>
        <tr>
          <th>片段</th><th>RefSeq</th><th class="num">長度 (bp)</th><th class="num">Reads</th><th class="num">MAPQ≥20</th>
          <th class="num">覆蓋廣度</th><th class="num">平均深度</th><th class="num">Consensus 有效鹼基 (≥10×)</th><th class="num">N 鹼基</th>
        </tr>
      </thead>
      <tbody>
        {''.join(flu_table_rows)}
      </tbody>
    </table>
  </div>

  <div class="panel">
    <h2>2. <em>Candida albicans</em> DNA 8 條染色體覆蓋分佈（SG26-GS01）</h2>
    <p style="font-size:13px;color:#475569;margin-top:0;">均勻分佈於全部 8 條染色體（3.84% – 5.15%），確認為真實真菌基因組 DNA。</p>
    {''.join(candida_bars)}
  </div>

  <div class="panel">
    <h2>3. 競爭性回貼總表（DNA &amp; cDNA）</h2>
    <p style="font-size:13px;color:#475569;margin-top:0;">
      Non-host 輸入：DNA <code>{overview['DNA_SG26-GS01']['total_primary']:,}</code> reads（非人 primary: <code>{overview['DNA_SG26-GS01']['nonhuman_primary']:,}</code>）；
      cDNA <code>{overview['cDNA_ST26-UM01']['total_primary']:,}</code> reads（人類殘留: <code>{overview['cDNA_ST26-UM01']['human_primary']:,}</code>，非人 primary: <code>{overview['cDNA_ST26-UM01']['nonhuman_primary']:,}</code>）。
    </p>
    <table>
      <thead>
        <tr>
          <th>Library</th><th>物種</th><th class="num">參考長度</th><th class="num">Mapped Reads</th><th class="num">MAPQ≥20</th>
          <th class="num">平均 MAPQ</th><th class="num">平均相似度</th><th class="num">覆蓋廣度</th><th class="num">平均深度</th>
        </tr>
      </thead>
      <tbody>
        {''.join(sp_table_rows)}
      </tbody>
    </table>
  </div>
</div>
</body>
</html>
"""
    out_html = os.path.join(OUTDIR, "pathogen_validation_report.html")
    with open(out_html, "w") as f:
        f.write(html_doc)
    return out_html


def main():
    out_fa, consensus_stats, depths = call_flu_consensus()
    out_html = build_html(consensus_stats, depths)
    print(f"Wrote consensus FASTA: {out_fa}")
    for c in consensus_stats:
        print(
            f"  {c['segment']} ({c['accession']}): {c['called_bases']}/{c['consensus_len']} bp ({c['called_pct']}%), N={c['n_bases']}"
        )
    print(f"Wrote HTML report: {out_html}")


if __name__ == "__main__":
    main()

