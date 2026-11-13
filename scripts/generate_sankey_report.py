#!/usr/bin/env python3
"""Generate standalone Sankey diagram (HTML) and final read disposition table (TSV)."""

import base64
import csv
import html
import json
import os

OUTDIR = "analysis/pathogen_validation"

# Exact mutually-exclusive read disposition accounting for all 30,833,825 raw reads
DISPOSITION_ROWS = [
    # DNA (SG26-GS01) = 17,989,222
    ("DNA_SG26-GS01", "QC_Failed", "chopper (Q<10 or <200bp)", "Filtered out", 165147),
    ("DNA_SG26-GS01", "Host_Removed", "hostile (human-masked map-ont)", "Host (Homo sapiens)", 17590626),
    ("DNA_SG26-GS01", "NonHost_FLU_H3N2", "Competitive minimap2", "Confirmed Pathogen (Influenza A H3N2)", 46),
    ("DNA_SG26-GS01", "NonHost_Candida_albicans", "Competitive minimap2", "Confirmed Fungal Signal (Candida albicans)", 332),
    ("DNA_SG26-GS01", "NonHost_Other_5Refs", "Competitive minimap2", "Low-MAPQ / Background (M. restricta, E. coli, etc.)", 117),
    ("DNA_SG26-GS01", "NonHost_Unmapped_K2_Classified", "Kraken2 PlusPF (unmapped to 8 refs)", "Residual Human / Bacterial k-mer / Kitome", 217970),
    ("DNA_SG26-GS01", "NonHost_Unmapped_K2_Unclassified", "Kraken2 + DIAMOND blastx (viral_refseq_prot)", "Screened Unknown (No Novel Viral Pathogen)", 14984),
    # cDNA (ST26-UM01) = 12,844,603
    ("cDNA_ST26-UM01", "QC_Failed", "chopper (Q<10 or <200bp)", "Filtered out", 1109174),
    ("cDNA_ST26-UM01", "Host_Removed", "hostile (human-masked splice)", "Host (Homo sapiens)", 4806406),
    ("cDNA_ST26-UM01", "NonHost_Residual_Human", "Competitive minimap2", "Residual Host (Homo sapiens)", 820534),
    ("cDNA_ST26-UM01", "NonHost_FLU_H3N2", "Competitive minimap2", "Confirmed Primary Pathogen (Influenza A H3N2)", 3694),
    ("cDNA_ST26-UM01", "NonHost_Candida_albicans", "Competitive minimap2", "Confirmed Fungal Signal (Candida albicans)", 1089),
    ("cDNA_ST26-UM01", "NonHost_Other_5Refs", "Competitive minimap2", "Low-MAPQ / Background (M. restricta, E. coli, etc.)", 2259),
    ("cDNA_ST26-UM01", "NonHost_Unmapped_K2_Classified", "Kraken2 PlusPF (unmapped to 8 refs)", "Residual Human / Bacterial k-mer / Kitome", 5698299),
    ("cDNA_ST26-UM01", "NonHost_Unmapped_TSO_Primer", "Motif check (unmapped K2_U)", "cDNA Library Primer Concatemer (TSO/SSP)", 402962),
    ("cDNA_ST26-UM01", "NonHost_Unmapped_Other_K2U", "Motif check + DIAMOND blastx (viral_refseq_prot)", "Screened Unknown (No Novel Viral Pathogen)", 186),
]


def write_tsv():
    total_all = sum(r[4] for r in DISPOSITION_ROWS)
    assert total_all == 30833825, f"Expected 30833825, got {total_all}"
    lib_totals = {}
    for lib, _, _, _, cnt in DISPOSITION_ROWS:
        lib_totals[lib] = lib_totals.get(lib, 0) + cnt

    out_tsv = os.path.join(OUTDIR, "read_disposition_summary.tsv")
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["library", "final_category", "method", "interpretation", "reads", "pct_of_library", "pct_of_total_raw"])
        for lib, cat, method, interp, cnt in DISPOSITION_ROWS:
            w.writerow([
                lib,
                cat,
                method,
                interp,
                cnt,
                f"{cnt / lib_totals[lib] * 100:.4f}",
                f"{cnt / total_all * 100:.4f}",
            ])
    return out_tsv, total_all, lib_totals


def render_sankey_svg(nodes, links, width=1120, height=540, title=""):
    """Render a clean multi-column Sankey SVG with sqrt-assisted scaling so small pathogen flows remain visible."""
    import math

    cols = sorted({n["col"] for n in nodes.values()})
    col_x = {c: 24 + i * ((width - 360) / max(1, len(cols) - 1)) for i, c in enumerate(cols)}
    node_w = 18
    pad_top = 44
    pad_bot = 20
    avail_h = height - pad_top - pad_bot

    # Compute node values
    for nid, n in nodes.items():
        in_v = sum(l["val"] for l in links if l["target"] == nid)
        out_v = sum(l["val"] for l in links if l["source"] == nid)
        n["val"] = max(in_v, out_v, n.get("val", 0))

    # Position nodes per column using power=0.45 scaling so 3,694 reads is clearly visible alongside 17M reads
    def scale_v(v):
        return math.pow(max(v, 1), 0.45)

    for c in cols:
        cnodes = [nid for nid, n in nodes.items() if n["col"] == c]
        gap = 12
        total_scaled = sum(scale_v(nodes[nid]["val"]) for nid in cnodes)
        usable_h = avail_h - gap * max(0, len(cnodes) - 1)
        y = pad_top
        for nid in cnodes:
            nh = max(6.0, usable_h * (scale_v(nodes[nid]["val"]) / total_scaled))
            nodes[nid]["x"] = col_x[c]
            nodes[nid]["y"] = y
            nodes[nid]["h"] = nh
            nodes[nid]["out_y"] = y
            nodes[nid]["in_y"] = y
            y += nh + gap

    # Compute outgoing/incoming link offsets proportional to node height
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
        tip = f"{s['label']} -> {t['label']}: {l['val']:,} reads ({l['val']/30833825*100:.3f}%)"
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
        anchor = "start"
        svg_nodes.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{node_w}" height="{h:.1f}" rx="3" fill="{color}">'
            f'<title>{html.escape(lbl)}</title></rect>'
            f'<text x="{tx:.1f}" y="{y + h/2 + 4:.1f}" font-size="11.5" font-weight="600" fill="#0f172a" text-anchor="{anchor}">{html.escape(lbl)}</text>'
        )

    return f"""<svg viewBox="0 0 {width} {height}" style="width:100%;height:auto;background:#fff;border:1px solid #e2e8f0;border-radius:10px;">
  <text x="20" y="24" font-size="14" font-weight="700" fill="#0f172a">{html.escape(title)}</text>
  {''.join(svg_links)}
  {''.join(svg_nodes)}
</svg>"""


def build_sankey_html(total_all, lib_totals):
    # Full pipeline Sankey (Raw -> Library -> QC/Host -> Final Disposition)
    nodes = {
        "raw": {"col": 0, "label": "All Raw Reads", "color": "#1e293b"},
        "dna": {"col": 1, "label": "DNA SG26-GS01", "color": "#2563eb"},
        "cdna": {"col": 1, "label": "cDNA ST26-UM01", "color": "#7c3aed"},
        "qc_fail": {"col": 2, "label": "QC Failed (chopper)", "color": "#94a3b8"},
        "host": {"col": 2, "label": "Host Removed (hostile)", "color": "#0284c7"},
        "nonhost": {"col": 2, "label": "Non-Host Retained", "color": "#0f766e"},
        "flu": {"col": 3, "label": "Influenza A (H3N2)", "color": "#16a34a"},
        "candida": {"col": 3, "label": "Candida albicans", "color": "#15803d"},
        "res_human": {"col": 3, "label": "Residual Human (mapped)", "color": "#0ea5e9"},
        "other5": {"col": 3, "label": "Other 5 Refs (Low-MAPQ)", "color": "#f59e0b"},
        "k2_other": {"col": 3, "label": "Unmapped K2 Classified", "color": "#d97706"},
        "tso": {"col": 3, "label": "TSO/SSP Primer Concatemers", "color": "#ea580c"},
        "k2_u": {"col": 3, "label": "Screened Unknown (DIAMOND: No Novel Virus)", "color": "#64748b"},
    }
    links = [
        {"source": "raw", "target": "dna", "val": 17989222, "color": "#2563eb"},
        {"source": "raw", "target": "cdna", "val": 12844603, "color": "#7c3aed"},
        {"source": "dna", "target": "qc_fail", "val": 165147, "color": "#94a3b8"},
        {"source": "cdna", "target": "qc_fail", "val": 1109174, "color": "#94a3b8"},
        {"source": "dna", "target": "host", "val": 17590626, "color": "#0284c7"},
        {"source": "cdna", "target": "host", "val": 4806406, "color": "#0284c7"},
        {"source": "dna", "target": "nonhost", "val": 233449, "color": "#0f766e"},
        {"source": "cdna", "target": "nonhost", "val": 6929023, "color": "#0f766e"},
        {"source": "nonhost", "target": "flu", "val": 3740, "color": "#16a34a"},
        {"source": "nonhost", "target": "candida", "val": 1421, "color": "#15803d"},
        {"source": "nonhost", "target": "res_human", "val": 820534, "color": "#0ea5e9"},
        {"source": "nonhost", "target": "other5", "val": 2376, "color": "#f59e0b"},
        {"source": "nonhost", "target": "k2_other", "val": 5916269, "color": "#d97706"},
        {"source": "nonhost", "target": "tso", "val": 402962, "color": "#ea580c"},
        {"source": "nonhost", "target": "k2_u", "val": 15170, "color": "#64748b"},
    ]
    sankey_svg = render_sankey_svg(
        nodes,
        links,
        width=1160,
        height=560,
        title="Sample_001 全序列流向與最終判定 Sankey 圖（總計 30,833,825 Raw Reads；採冪次比例縮放以兼顧微量病原體可視性）",
    )

    table_rows = []
    for lib, cat, method, interp, cnt in DISPOSITION_ROWS:
        pct_lib = cnt / lib_totals[lib] * 100.0
        pct_all = cnt / total_all * 100.0
        hl = " style='background:#f0fdf4;font-weight:600;'" if "Confirmed" in interp else ""
        table_rows.append(
            f"<tr{hl}><td>{html.escape(lib)}</td><td><code>{html.escape(cat)}</code></td>"
            f"<td>{html.escape(method)}</td><td>{html.escape(interp)}</td>"
            f"<td class='num'>{cnt:,}</td><td class='num'>{pct_lib:.4f}%</td><td class='num'>{pct_all:.4f}%</td></tr>"
        )

    logo_b64 = ""
    logo_path = os.path.join("webapp", "assets", "clinmeta-logo-light.jpg")
    if os.path.exists(logo_path):
        with open(logo_path, "rb") as lf:
            logo_b64 = base64.b64encode(lf.read()).decode("ascii")
    logo_html = (
        f'<img src="data:image/jpeg;base64,{logo_b64}" alt="ClinMeta" style="height:56px;border-radius:6px;">'
        if logo_b64 else ""
    )

    out_html = os.path.join(OUTDIR, "read_disposition_sankey.html")
    doc = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<title>Sample_001 — 全序列最終判定 Sankey 圖與統計表</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; background: #f1f5f9; color: #0f172a; line-height: 1.5; }}
  .container {{ max-width: 1220px; margin: 28px auto; padding: 0 20px; }}
  h1 {{ font-size: 24px; margin-bottom: 4px; }}
  .sub {{ color: #475569; margin-bottom: 20px; font-size: 14px; }}
  .panel {{ background: #fff; border-radius: 10px; padding: 20px; margin-bottom: 22px; box-shadow: 0 1px 3px rgba(15,23,42,0.08); }}
  .panel h2 {{ font-size: 18px; margin-top: 0; margin-bottom: 12px; border-bottom: 1px solid #e2e8f0; padding-bottom: 8px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
  th, td {{ padding: 8px 10px; border-bottom: 1px solid #e2e8f0; text-align: left; }}
  th {{ background: #f8fafc; font-weight: 600; color: #334155; }}
  td.num, th.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  code {{ background: #f1f5f9; padding: 2px 5px; border-radius: 4px; font-size: 12px; }}
  .sankey-link:hover {{ stroke-opacity: 0.78; }}
</style>
</head>
<body>
<div class="container">
  <div style="display:flex;align-items:center;justify-content:space-between;gap:16px;flex-wrap:wrap;margin-bottom:14px;background:#fff;padding:14px 20px;border-radius:10px;box-shadow:0 1px 3px rgba(15,23,42,0.08);">
    <div>
      <h1 style="margin:0 0 4px;">ClinMeta <code>sample_001</code> — 全序列最終判定 Sankey 圖與明細表</h1>
      <div class="sub" style="margin:0;">完整追蹤全數 <strong>{total_all:,}</strong> 條原始定序 reads（DNA <code>SG26-GS01</code> {lib_totals['DNA_SG26-GS01']:,} + cDNA <code>ST26-UM01</code> {lib_totals['cDNA_ST26-UM01']:,}）從 QC、去宿主到最終病原體／背景／人為產物判定的互斥歸屬。</div>
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
    const SK_TOTAL = {total_all};
    function renderInteractiveSankey() {{
      const scope = document.getElementById('sk-scope').value;
      const scaleMode = document.getElementById('sk-scale').value;
      const lblMode = document.getElementById('sk-lbl').value;
      const nodes = JSON.parse(JSON.stringify(SK_RAW_NODES));
      let links = JSON.parse(JSON.stringify(SK_RAW_LINKS));
      if (scope === 'nonhost') {{
        const keep = new Set(['nonhost', 'flu', 'candida', 'res_human', 'other5', 'k2_other', 'tso', 'k2_u']);
        links = links.filter(l => keep.has(l.source) && keep.has(l.target));
        Object.keys(nodes).forEach(k => {{ if (!keep.has(k)) delete nodes[k]; else if (nodes[k].col >= 2) nodes[k].col -= 2; }});
      }}
      const W = 1160, H = 560, nodeW = 18, padTop = 46, padBot = 22, availH = H - padTop - padBot;
      const cols = [...new Set(Object.values(nodes).map(n => n.col))].sort((a,b) => a-b);
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
      const title = scope === 'nonhost' ? 'Sample_001 非人類 (Non-Host) 序列流向 Sankey 圖' : 'Sample_001 全序列流向與最終判定 Sankey 圖（總計 30,833,825 Raw Reads）';
      document.getElementById('sankey-box').innerHTML = `<svg viewBox="0 0 ${{W}} ${{H}}" style="width:100%;height:auto;background:#fff;border:1px solid #e2e8f0;border-radius:10px;"><text x="20" y="26" font-size="14" font-weight="700" fill="#0f172a">${{esc(title)}}</text>${{paths}}${{rects}}</svg>`;
    }}
    function downloadSankeySvg() {{
      const svg = document.querySelector('#sankey-box svg');
      if (!svg) return;
      const blob = new Blob([svg.outerHTML], {{type:'image/svg+xml;charset=utf-8'}});
      const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'sample_001_sankey.svg'; a.click();
    }}
    renderInteractiveSankey();
  </script>

  <div class="panel">
    <h2>全序列互斥最終判定明細表（合計 {total_all:,} reads = 100.0000%）</h2>
    <table>
      <thead>
        <tr>
          <th>Library</th><th>最終分類項目</th><th>判定方法</th><th>生物／臨床意義</th>
          <th class="num">Reads</th><th class="num">佔該 Library (%)</th><th class="num">佔總 Raw (%)</th>
        </tr>
      </thead>
      <tbody>
        {''.join(table_rows)}
      </tbody>
    </table>
  </div>
</div>
</body>
</html>
"""
    with open(out_html, "w") as f:
        f.write(doc)
    return out_html


def main():
    out_tsv, total_all, lib_totals = write_tsv()
    out_html = build_sankey_html(total_all, lib_totals)
    print(f"Wrote {out_tsv} (total={total_all:,})")
    print(f"Wrote {out_html}")


if __name__ == "__main__":
    main()

