#!/usr/bin/env python3
"""Build a single-file interactive HTML report (plotly.js, incl. Sankey) from Kraken2 reports.

Usage: python scripts/kraken2_report_html.py [clean_results/kraken2] [output.html]
The plotly.js bundle (scripts/vendor/plotly-4.3.0.min.js, MIT) is inlined so the report
opens offline by double-click; no server or internet access is needed.
"""
import csv
import datetime as dt
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from summarize_kraken2 import CONTAMINANT_GENERA, MIN_READS_REPLICATE  # noqa: E402

PLOTLY = HERE / 'vendor' / 'plotly-4.3.0.min.js'


def parse_report(path):
    nodes, stack = [], []
    for line in path.read_text().splitlines():
        f = line.split('\t')
        if len(f) < 8:
            continue
        raw = f[7]
        depth = (len(raw) - len(raw.lstrip(' '))) // 2
        rank = f[5]
        node = dict(t=f[6], n=raw.strip(), r=rank, c=int(f[1]), d=int(f[2]),
                    m=int(f[4]), p=-1)
        if rank == 'U':
            nodes.append(node)
            continue
        while stack and stack[-1][0] >= depth:
            stack.pop()
        node['p'] = stack[-1][1] if stack else -1
        nodes.append(node)
        stack.append((depth, len(nodes) - 1))
    return nodes


def read_tsv(path):
    if not path.is_file():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def main(kdir, html):
    kdir = Path(kdir).resolve()
    results = kdir.parent
    sources = read_tsv(results / 'source_summary.tsv')
    src_info = {f"source.{r['source_id']}": r for r in sources}
    sample = {}
    if (results / 'sample_summary.json').is_file():
        sample = json.loads((results / 'sample_summary.json').read_text())
    reports = []
    for rep in sorted(kdir.glob('*.report.txt')):
        name = rep.name[:-len('.report.txt')]
        nodes = parse_report(rep)
        total = sum(n['c'] for n in nodes if n['r'] in ('U', 'R'))
        unclassified = sum(n['c'] for n in nodes if n['r'] == 'U')
        if name in src_info:
            s = src_info[name]
            run_id = s.get('run_id') or s.get('source_id') or ''
            lib_type, label = s.get('library_type', 'Unknown'), f"{s.get('library_type', '')} · {run_id}".strip(' ·')
            kind = 'source'
        else:
            lib = name.split('.')[-1]
            lib_type = next((r['library_type'] for r in sources if r['library_id'] == lib), lib)
            label, kind = f'{lib_type} 合併 ({lib})', 'library'
        reports.append(dict(name=name, label=label, kind=kind, lib_type=lib_type,
                            total=total, unclassified=unclassified, nodes=nodes))
    if not reports:
        raise SystemExit(f'No *.report.txt in {kdir}')
    # Species-level flags: contaminant genus list + presence in both sequencing runs.
    flags = {}
    for rep in reports:
        for n in rep['nodes']:
            if n['r'] == 'S':
                flags.setdefault(n['t'], dict(name=n['n'], runs={}))
    for taxid, f in flags.items():
        f['contam'] = f['name'].split()[0] in CONTAMINANT_GENERA
        for rep in reports:
            if rep['kind'] != 'source':
                continue
            reads = next((n['c'] for n in rep['nodes'] if n['t'] == taxid and n['r'] == 'S'), 0)
            f['runs'].setdefault(rep['lib_type'], []).append(reads)
        f['both'] = {t: all(x >= MIN_READS_REPLICATE for x in v) for t, v in f['runs'].items()}
        del f['runs']
    meta = dict(
        generated=dt.datetime.now().astimezone().isoformat(timespec='seconds'),
        database=(kdir / 'database_path.txt').read_text().strip()
        if (kdir / 'database_path.txt').is_file() else 'unknown',
        kraken2=(kdir / 'kraken2_version.txt').read_text().strip()
        if (kdir / 'kraken2_version.txt').is_file() else 'unknown',
        sample_id=sample.get('sample_id', ''), mode=sample.get('mode', ''),
        config={k: sample.get('config', {}).get(k) for k in
                ('quality', 'min_length', 'headcrop', 'tailcrop')},
        min_reads_replicate=MIN_READS_REPLICATE)
    data = dict(meta=meta, sources=sources, reports=reports,
                flags={k: dict(contam=v['contam'], both=v['both']) for k, v in flags.items()})
    payload = json.dumps(data, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    page = TEMPLATE.replace('__PLOTLY__', PLOTLY.read_text()).replace('__DATA__', payload)
    out = Path(html) if html else kdir / 'report.html'
    out.write_text(page, encoding='utf-8')
    print(f'Wrote {out} ({out.stat().st_size / 1e6:.1f} MB)')


TEMPLATE = r'''<!DOCTYPE html>
<html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Nanopore 去宿主後微生物分類報告</title>
<style>
:root{--fg:#1f2937;--muted:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--card:#fff;--accent:#2563eb}
*{box-sizing:border-box}body{margin:0;font-family:system-ui,-apple-system,"Noto Sans TC","Microsoft JhengHei",sans-serif;color:var(--fg);background:var(--bg)}
header{background:#0f172a;color:#fff;padding:20px 32px}header h1{margin:0 0 4px;font-size:22px}header p{margin:0;color:#cbd5e1;font-size:13px}
nav{position:sticky;top:0;z-index:5;background:#fff;border-bottom:1px solid var(--line);padding:8px 32px;display:flex;gap:18px;flex-wrap:wrap}
nav a{color:var(--accent);text-decoration:none;font-size:14px}
main{max-width:1400px;margin:0 auto;padding:20px 32px 60px}
section{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px 20px;margin-bottom:20px}
h2{font-size:18px;margin:0 0 10px}.muted{color:var(--muted);font-size:13px}
.warn{background:#fff7ed;border-left:4px solid #f97316;padding:10px 14px;border-radius:6px;font-size:14px;line-height:1.6}
.controls{display:flex;gap:16px;flex-wrap:wrap;align-items:center;margin:8px 0 12px;font-size:14px}
.controls label{display:flex;gap:6px;align-items:center}
select,input[type=number],input[type=search]{padding:4px 6px;border:1px solid #cbd5e1;border-radius:6px;font-size:14px}
table{border-collapse:collapse;width:100%;font-size:13px}th,td{border-bottom:1px solid var(--line);padding:6px 8px;text-align:right;white-space:nowrap}
th{background:#f1f5f9;cursor:pointer;position:sticky;top:0}td:first-child,th:first-child,td.l,th.l{text-align:left}
.tablewrap{max-height:560px;overflow:auto;border:1px solid var(--line);border-radius:6px}
.tag{display:inline-block;padding:1px 6px;border-radius:9px;font-size:11px;margin-left:4px}
.tag.c{background:#fee2e2;color:#991b1b}.tag.ok{background:#dcfce7;color:#166534}.tag.no{background:#f1f5f9;color:#64748b}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px}@media(max-width:1000px){.grid2{grid-template-columns:1fr}}
</style>
<script>__PLOTLY__</script>
</head><body>
<header><h1>Nanopore 去宿主後微生物分類報告</h1><p id="meta"></p></header>
<nav><a href="#s-qc">QC／去宿主</a><a href="#s-cls">分類概況</a><a href="#s-sankey">Sankey</a><a href="#s-heat">物種熱圖</a><a href="#s-bar">Top 物種</a><a href="#s-table">物種表</a><a href="#s-notes">判讀說明</a></nav>
<main>
<section><div class="warn"><b>判讀限制：</b>本批次<b>沒有陰性對照</b>；Kraken2 結果為 k-mer 分類的<b>候選微生物證據</b>，不等於感染診斷。
低量訊號無法排除試劑污染或 barcode 串擾；<span class="tag c">常見污染屬</span>標記代表需謹慎判讀，並非排除（多個屬同時為病原菌）。
候選物種需再以參考基因組回貼、覆蓋廣度、與人類／近緣物種競爭比對確認。DNA 與 cDNA 的百分比不可直接比較。</div></section>

<section id="s-qc"><h2>QC 與去宿主（逐批來源）</h2><p class="muted">raw → chopper QC → hostile 去宿主後非人類 reads。y 軸為對數刻度。</p>
<div class="grid2"><div id="qcBar" style="height:420px"></div><div id="qcTable"></div></div></section>

<section id="s-cls"><h2>分類概況</h2><p class="muted">各 Kraken2 報告中已分類、未分類，以及仍被分類為 <i>Homo sapiens</i> 的殘留 reads 比例（分母＝該報告的非人類 reads）。</p>
<div id="clsBar" style="height:420px"></div></section>

<section id="s-sankey"><h2>Sankey：分類階層流向</h2>
<div class="controls">
<label>報告 <select id="skReport"></select></label>
<label>Top 物種數 <input id="skTop" type="number" value="25" min="5" max="200" step="5" style="width:70px"></label>
<label>物種最少 reads <input id="skMin" type="number" value="3" min="1" style="width:70px"></label>
<label><input id="skHuman" type="checkbox" checked>排除 Homo sapiens</label>
<label><input id="skContam" type="checkbox">排除常見污染屬</label>
<span>階層：</span><span id="skRanks"></span></div>
<p class="muted">流量＝所選 Top 物種的 clade reads 往上加總；節點滑鼠停留可看該分類群在報告中的實際 clade reads。顏色依 domain：細菌藍、病毒紅、真核綠、古菌紫。</p>
<div id="sankey" style="height:760px"></div></section>

<section id="s-heat"><h2>物種熱圖（跨報告）</h2>
<div class="controls"><label>物種數 <input id="hmTop" type="number" value="40" min="10" max="200" step="10" style="width:70px"></label>
<label><input id="hmHuman" type="checkbox" checked>排除 Homo sapiens</label>
<label>色階 <select id="hmScale"><option value="log">log10(% of non-host)</option><option value="lin">% of non-host</option></select></label></div>
<div id="heatmap" style="height:900px"></div></section>

<section id="s-bar"><h2>Top 物種（單一報告）</h2>
<div class="controls"><label>報告 <select id="barReport"></select></label>
<label>物種數 <input id="barTop" type="number" value="30" min="5" max="100" step="5" style="width:70px"></label>
<label><input id="barHuman" type="checkbox" checked>排除 Homo sapiens</label></div>
<p class="muted">紅色＝常見試劑污染屬；顏色深淺不代表可信度。hover 顯示 distinct minimizers（越多代表 reads 分布於基因組更多區段）。</p>
<div id="topBar" style="height:760px"></div></section>

<section id="s-table"><h2>物種表</h2>
<div class="controls"><label>搜尋 <input id="tblSearch" type="search" placeholder="物種名稱 / taxid"></label>
<label>最少 reads（任一報告）<input id="tblMin" type="number" value="1" min="0" style="width:80px"></label>
<label><input id="tblBoth" type="checkbox">只顯示兩次定序皆偵測</label>
<label><input id="tblNoContam" type="checkbox">隱藏常見污染屬</label>
<span class="muted" id="tblCount"></span></div>
<div class="tablewrap"><table id="tbl"></table></div></section>

<section id="s-notes"><h2>判讀說明與參數</h2><div id="notes" class="muted" style="line-height:1.7"></div></section>
</main>
<script id="data" type="application/json">__DATA__</script>
<script>
const D = JSON.parse(document.getElementById('data').textContent);
const R = D.reports, F = D.flags, M = D.meta;
const RANK_LABEL = {D:'Domain',K:'Kingdom',P:'Phylum',C:'Class',O:'Order',F:'Family',G:'Genus',S:'Species'};
const DOM_COLOR = {Bacteria:'#2563eb',Viruses:'#dc2626',Eukaryota:'#16a34a',Archaea:'#7c3aed'};
const fmt = n => n.toLocaleString('en-US');
const pct = (a,b) => b ? (100*a/b) : 0;
const cfg = {responsive:true, displaylogo:false, toImageButtonOptions:{format:'svg'}};
const isHuman = n => n.n === 'Homo sapiens';
const isContam = n => (F[n.t] && F[n.t].contam);

document.getElementById('meta').textContent =
  `樣本 ${M.sample_id||'—'} · 模式 ${M.mode||'—'} · ${M.kraken2} · DB: ${M.database.split('/').pop()} · 產生時間 ${M.generated}`;

// helpers ---------------------------------------------------------------
function domainOf(rep, i){ let n=rep.nodes[i]; while(n){ if(n.r==='D') return n.n; if(n.p<0) return null; n=rep.nodes[n.p]; } return null; }
function speciesOf(rep){ return rep.nodes.map((n,i)=>({...n,i})).filter(n=>n.r==='S'); }
function fillSelect(id){ const s=document.getElementById(id); R.forEach((r,i)=>{const o=document.createElement('option');o.value=i;o.textContent=r.label;s.appendChild(o);});
  const lib=R.findIndex(r=>r.kind==='library'); if(lib>=0) s.value=lib; }
function tags(t){ const f=F[t]||{}; let h=f.contam?'<span class="tag c">常見污染屬</span>':'';
  for(const [k,v] of Object.entries(f.both||{})) h+=`<span class="tag ${v?'ok':'no'}">${k} 兩次${v?'皆有':'未皆有'}</span>`; return h; }

// QC --------------------------------------------------------------------
(function(){
  const S=D.sources; if(!S.length){document.getElementById('qcBar').innerHTML='<p class="muted">無 source_summary.tsv</p>';return;}
  const x=S.map(s=>`${s.library_type} · ${s.run_id}`);
  const tr=(k,name,c)=>({type:'bar',name,x,y:S.map(s=>+s[k]),marker:{color:c},hovertemplate:'%{x}<br>'+name+': %{y:,}<extra></extra>'});
  Plotly.newPlot('qcBar',[tr('raw_reads','Raw','#94a3b8'),tr('qc_reads','QC 後','#60a5fa'),tr('nonhost_reads','非人類','#f97316')],
    {barmode:'group',yaxis:{type:'log',title:'reads'},margin:{t:24,l:70,r:160,b:60},legend:{orientation:'v',x:1.02,y:1}},cfg);
  let h='<table><tr><th class="l">來源</th><th>Raw reads</th><th>Raw bases</th><th>QC 後</th><th>移除人類</th><th>非人類 reads</th><th>非人類 %</th><th>非人類 bases</th></tr>';
  S.forEach((s,i)=>{h+=`<tr><td>${x[i]}<br><span class="muted">${s.library_id}</span></td><td>${fmt(+s.raw_reads)}</td><td>${fmt(+s.raw_bases)}</td><td>${fmt(+s.qc_reads)}</td><td>${fmt(+s.host_removed_reads)}</td><td>${fmt(+s.nonhost_reads)}</td><td>${pct(+s.nonhost_reads,+s.qc_reads).toFixed(2)}%</td><td>${fmt(+s.nonhost_bases)}</td></tr>`;});
  document.getElementById('qcTable').innerHTML=h+'</table>';
})();

// Classification overview ----------------------------------------------
(function(){
  const x=R.map(r=>r.label);
  const human=R.map(r=>{const n=r.nodes.find(isHuman);return n?n.c:0;});
  const cls=R.map((r,i)=>r.total-r.unclassified-human[i]);
  const mk=(y,name,c,abs)=>({type:'bar',name,x,y:y.map((v,i)=>pct(v,R[i].total)),customdata:y,marker:{color:c},
    hovertemplate:'%{x}<br>'+name+': %{customdata:,} reads (%{y:.2f}%)<extra></extra>'});
  Plotly.newPlot('clsBar',[mk(cls,'已分類（非人類）','#2563eb'),mk(human,'殘留 Homo sapiens','#a3a3a3'),mk(R.map(r=>r.unclassified),'未分類','#e5e7eb')],
    {barmode:'stack',yaxis:{title:'% of non-host reads',range:[0,100]},margin:{t:24,l:60,r:170,b:80},legend:{orientation:'v',x:1.02,y:1}},cfg);
})();

// Sankey ----------------------------------------------------------------
const ALL_RANKS=['D','P','C','O','F','G','S'], DEFAULT_RANKS=new Set(['D','P','F','G','S']);
(function(){ const box=document.getElementById('skRanks');
  ALL_RANKS.forEach(r=>{const l=document.createElement('label');l.innerHTML=`<input type="checkbox" value="${r}" ${DEFAULT_RANKS.has(r)?'checked':''} ${r==='S'?'disabled':''}>${RANK_LABEL[r]}`;box.appendChild(l);}); })();
function drawSankey(){
  const rep=R[+document.getElementById('skReport').value];
  const top=+document.getElementById('skTop').value, min=+document.getElementById('skMin').value;
  const exH=document.getElementById('skHuman').checked, exC=document.getElementById('skContam').checked;
  const ranks=new Set([...document.querySelectorAll('#skRanks input:checked')].map(e=>e.value)); ranks.add('S');
  let sp=speciesOf(rep).filter(n=>n.c>=min && !(exH&&isHuman(n)) && !(exC&&isContam(n))).sort((a,b)=>b.c-a.c).slice(0,top);
  const idx=new Map(), labels=[], colors=[], custom=[]; const links=new Map();
  const nodeId=(i,label,dom,clade,rank)=>{ if(!idx.has(i)){ idx.set(i,labels.length); labels.push(label);
      colors.push(DOM_COLOR[dom]||'#64748b'); custom.push([fmt(clade),RANK_LABEL[rank]||rank]); } return idx.get(i); };
  const rootLabel='非人類 reads（'+rep.label+'）';
  sp.forEach(s=>{
    const dom=domainOf(rep,s.i); const chain=[]; let j=s.i;
    while(j>=0){ const n=rep.nodes[j]; if(ranks.has(n.r)) chain.unshift(j); j=n.p; }
    let prev=nodeId(-1,rootLabel,null,rep.total,'root');
    chain.forEach(k=>{ const n=rep.nodes[k]; const cur=nodeId(k,n.n,dom,n.c,n.r);
      const key=prev+'>'+cur; links.set(key,(links.get(key)||0)+s.c); prev=cur; });
  });
  const L=[...links.entries()].map(([k,v])=>{const [a,b]=k.split('>').map(Number);return {a,b,v};});
  const hex2rgba=(h,a)=>{const n=parseInt(h.slice(1),16);return `rgba(${n>>16&255},${n>>8&255},${n&255},${a})`;};
  Plotly.react('sankey',[{type:'sankey',arrangement:'snap',valueformat:',d',
    node:{label:labels,color:colors,customdata:custom,pad:10,thickness:16,line:{width:0},
      hovertemplate:'%{label}<br>%{customdata[1]}<br>clade reads（報告）: %{customdata[0]}<br>Sankey 流量: %{value:,}<extra></extra>'},
    link:{source:L.map(l=>l.a),target:L.map(l=>l.b),value:L.map(l=>l.v),color:L.map(l=>hex2rgba(colors[l.b],0.35)),
      hovertemplate:'%{source.label} → %{target.label}<br>%{value:,} reads<extra></extra>'}}],
    {margin:{t:10,l:10,r:10,b:10},font:{size:11},
     annotations: sp.length?[]:[{text:'沒有符合條件的物種',showarrow:false,x:0.5,y:0.5,xref:'paper',yref:'paper'}]},cfg);
}

// Heatmap ---------------------------------------------------------------
function drawHeat(){
  const top=+document.getElementById('hmTop').value, exH=document.getElementById('hmHuman').checked, lin=document.getElementById('hmScale').value==='lin';
  const best=new Map();
  R.forEach(r=>speciesOf(r).forEach(n=>{ if(exH&&isHuman(n)) return; const p=pct(n.c,r.total);
    const b=best.get(n.t)||{name:n.n,max:0}; b.max=Math.max(b.max,p); best.set(n.t,b); }));
  const taxa=[...best.entries()].sort((a,b)=>b[1].max-a[1].max).slice(0,top).reverse();
  const z=[],text=[];
  taxa.forEach(([t])=>{ const row=[],tr=[]; R.forEach(r=>{ const n=r.nodes.find(x=>x.t===t&&x.r==='S'); const c=n?n.c:0, p=pct(c,r.total);
      row.push(c? (lin?p:Math.log10(p)) : null); tr.push(`${fmt(c)} reads · ${p.toFixed(3)}%<br>distinct minimizers: ${n?fmt(n.m):0}`); }); z.push(row); text.push(tr); });
  const y=taxa.map(([t,b])=>(F[t]&&F[t].contam?'⚠ ':'')+b.name);
  Plotly.react('heatmap',[{type:'heatmap',z,x:R.map(r=>r.label),y,text,colorscale:'YlOrRd',reversescale:false,hoverongaps:false,
    colorbar:{title:{text:lin?'%':'log10(%)'}},hovertemplate:'%{y}<br>%{x}<br>%{text}<extra></extra>'}],
    {margin:{t:10,l:280,r:10,b:110},xaxis:{tickangle:-25},yaxis:{automargin:true,tickfont:{size:11}}},cfg);
}

// Top bar ---------------------------------------------------------------
function drawBar(){
  const rep=R[+document.getElementById('barReport').value], top=+document.getElementById('barTop').value, exH=document.getElementById('barHuman').checked;
  const sp=speciesOf(rep).filter(n=>!(exH&&isHuman(n))).sort((a,b)=>b.c-a.c).slice(0,top).reverse();
  Plotly.react('topBar',[{type:'bar',orientation:'h',y:sp.map(n=>n.n),x:sp.map(n=>n.c),
    marker:{color:sp.map(n=>isContam(n)?'#ef4444':(DOM_COLOR[domainOf(rep,n.i)]||'#64748b'))},
    customdata:sp.map(n=>[pct(n.c,rep.total).toFixed(3),fmt(n.m),fmt(n.d)]),
    hovertemplate:'%{y}<br>clade reads: %{x:,} (%{customdata[0]}%)<br>direct reads: %{customdata[2]}<br>distinct minimizers: %{customdata[1]}<extra></extra>'}],
    {margin:{t:10,l:260,r:20,b:50},xaxis:{title:'clade reads',type:'log'},yaxis:{automargin:true,tickfont:{size:11}}},cfg);
}

// Table -----------------------------------------------------------------
let sortKey=null, sortDir=-1;
function drawTable(){
  const q=document.getElementById('tblSearch').value.trim().toLowerCase(), min=+document.getElementById('tblMin').value;
  const both=document.getElementById('tblBoth').checked, noC=document.getElementById('tblNoContam').checked;
  const rows=new Map();
  R.forEach((r,ri)=>speciesOf(r).forEach(n=>{ const o=rows.get(n.t)||{t:n.t,name:n.n,v:R.map(()=>[0,0,0])}; o.v[ri]=[n.c,pct(n.c,r.total),n.m]; rows.set(n.t,o); }));
  let list=[...rows.values()].filter(o=>{ const f=F[o.t]||{};
    if(q && !o.name.toLowerCase().includes(q) && !o.t.includes(q)) return false;
    if(Math.max(...o.v.map(v=>v[0]))<min) return false;
    if(both && !Object.values(f.both||{}).some(Boolean)) return false;
    if(noC && f.contam) return false; return true; });
  const key=sortKey??R.findIndex(r=>r.kind==='library');
  list.sort((a,b)=> key==='name'? sortDir*a.name.localeCompare(b.name) : sortDir*(a.v[Math.max(key,0)][0]-b.v[Math.max(key,0)][0]));
  const shown=list.slice(0,500);
  let h='<tr><th class="l" data-k="name">物種</th><th class="l">taxid</th><th class="l">標記</th>'+R.map((r,i)=>`<th data-k="${i}">${r.label}<br><span class="muted">reads / % / distinct mm</span></th>`).join('')+'</tr>';
  shown.forEach(o=>{ h+=`<tr><td><i>${o.name}</i></td><td class="l">${o.t}</td><td class="l">${tags(o.t)}</td>`+o.v.map(v=>`<td>${v[0]?fmt(v[0])+' / '+v[1].toFixed(3)+'% / '+fmt(v[2]):'—'}</td>`).join('')+'</tr>'; });
  const t=document.getElementById('tbl'); t.innerHTML=h;
  document.getElementById('tblCount').textContent=`${fmt(list.length)} 個物種${list.length>500?'（顯示前 500）':''}`;
  t.querySelectorAll('th[data-k]').forEach(th=>th.onclick=()=>{ const k=th.dataset.k==='name'?'name':+th.dataset.k; sortDir=(sortKey===k)?-sortDir:-1; sortKey=k; drawTable(); });
}

// Notes -----------------------------------------------------------------
document.getElementById('notes').innerHTML = `
<ul>
<li>Kraken2 資料庫：<code>${M.database}</code>；${M.kraken2}。百分比分母＝該報告輸入的非人類 reads 總數（含未分類）。</li>
<li>QC 參數：Q≥${M.config.quality??'—'}、長度≥${M.config.min_length??'—'} bp、頭尾裁切 ${M.config.headcrop??'—'}/${M.config.tailcrop??'—'} bp。</li>
<li>「兩次皆有」＝同一 library 的兩次定序（技術重複）各至少 ${M.min_reads_replicate} reads。</li>
<li>「常見污染屬」依 Salter et al. 2014 (BMC Biol 12:87) 與 Eisenhofer et al. 2019 (Trends Microbiol 27:105)；僅為警示。</li>
<li>distinct minimizers 低而 reads 高，常代表 reads 集中於少數保守區域或低複雜度序列，需回貼確認。</li>
<li>Sankey 流量只計入所選 Top 物種；未達物種層級的 reads 不在圖中，請參考分類概況與物種表。</li>
<li>圖表右上角相機圖示可匯出 SVG。</li></ul>`;

// wire up ---------------------------------------------------------------
['skReport','barReport'].forEach(fillSelect);
['skReport','skTop','skMin','skHuman','skContam'].forEach(id=>document.getElementById(id).addEventListener('change',drawSankey));
document.getElementById('skRanks').addEventListener('change',drawSankey);
['hmTop','hmHuman','hmScale'].forEach(id=>document.getElementById(id).addEventListener('change',drawHeat));
['barReport','barTop','barHuman'].forEach(id=>document.getElementById(id).addEventListener('change',drawBar));
['tblSearch','tblMin','tblBoth','tblNoContam'].forEach(id=>document.getElementById(id).addEventListener('input',drawTable));
drawSankey(); drawHeat(); drawBar(); drawTable();
window.addEventListener('resize',()=>{['qcBar','clsBar','sankey','heat','barPlot'].forEach(id=>{const el=document.getElementById(id);if(el&&el.data)try{Plotly.Plots.resize(el);}catch(e){}});});
</script></body></html>
'''

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'clean_results/kraken2',
         sys.argv[2] if len(sys.argv) > 2 else None)

