#!/usr/bin/env python3
"""Audit every dehosted read using NCBI DUST and transparent complexity flags."""
import concurrent.futures,csv,gzip,json,math,re,subprocess
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/pathogen_validation'
DB=ROOT/'databases/kraken2/k2_pluspf_20261022'
parents={};ranks={}
for line in (DB/'nodes.dmp').open():
 f=line.split('|');parents[int(f[0])]=int(f[1]);ranks[int(f[0])]=f[2].strip()
species_cache={}
def species(t):
 if t in species_cache:return species_cache[t]
 orig=t;seen=set()
 while t not in seen:
  if ranks.get(t)=='species':species_cache[orig]=t;return t
  seen.add(t);t=parents.get(t,1)
 species_cache[orig]=0;return 0
rows=list(csv.DictReader((ROOT/'clean_results/source_summary.tsv').open(),delimiter='\t'))
def process(row):
 sid=row['source_id'];work=OUT/'sources'/sid;work.mkdir(parents=True,exist_ok=True)
 done=work/'complexity_summary.json'
 if done.exists():return json.loads(done.read_text())
 classes={}
 kraken=ROOT/'clean_results/kraken2'/f'source.{sid}.reads.kraken.gz'
 with gzip.open(kraken,'rt') as handle:
  for line in handle:
   f=line.split('\t',4)
   if f[0]=='C':
    tax=int(re.search(r'\(taxid (\d+)\)',f[2]).group(1))
    classes[f[1]]=(tax,species(tax))
   else:classes[f[1]]=(0,0)
 counts=Counter();taxstats=defaultdict(Counter)
 fields=['read_id','length','dust_bases','dust_fraction','entropy','dominant_fraction','low_complexity','kraken_taxid','kraken_species_taxid']
 with (work/'dustmasker.log').open('w') as log, gzip.open(work/'read_features.tsv.gz','wt') as features, \
      gzip.open(work/'informative.fa.gz','wt') as informative:
  writer=csv.writer(features,delimiter='\t');writer.writerow(fields)
  first=subprocess.Popen(['seqkit','fq2fa','-w','0',str(ROOT/'clean_results/sources'/sid/'clean.fastq.gz')],stdout=subprocess.PIPE,stderr=log)
  second=subprocess.Popen(['dustmasker','-outfmt','fasta'],stdin=first.stdout,stdout=subprocess.PIPE,stderr=log,text=True)
  first.stdout.close()
  def consume(header,parts):
   if header is None:return
   rid=header.split()[0];seq=''.join(parts);n=len(seq);upper=seq.upper()
   # DUST soft-mask fraction; global entropy and base dominance supplement local masking.
   masked=sum(c.islower() for c in seq);base_counts=[upper.count(c) for c in 'ACGT']
   entropy=-sum((c/n)*math.log2(c/n) for c in base_counts if c)
   dominant=max(base_counts)/n if n else 1
   low=(masked/n>=0.50 or entropy<1.2 or dominant>=0.8 or n-masked<100)
   tax,sp=classes.pop(rid)
   writer.writerow([rid,n,masked,f'{masked/n:.6f}',f'{entropy:.5f}',f'{dominant:.5f}',int(low),tax,sp])
   counts['reads']+=1;counts['bases']+=n;counts['dust_bases']+=masked
   counts['low_complexity_reads']+=low;counts['informative_reads']+=not low
   taxstats[sp]['reads']+=1;taxstats[sp]['low_complexity_reads']+=low
   taxstats[sp]['bases']+=n;taxstats[sp]['dust_bases']+=masked
   if not low:informative.write('>'+rid+'\n'+upper+'\n')
  header=None;parts=[]
  for line in second.stdout:
   if line.startswith('>'):
    consume(header,parts);header=line[1:].strip();parts=[]
   else:parts.append(line.strip())
  consume(header,parts)
  if second.wait()!=0 or first.wait()!=0:raise RuntimeError('DUST failed: '+sid)
 if classes:raise RuntimeError('Kraken/FASTA ID mismatch: '+sid)
 if counts['reads']!=int(row['nonhost_reads']):raise RuntimeError('Read count mismatch: '+sid)
 with (work/'complexity_by_species.tsv').open('w') as handle:
  w=csv.writer(handle,delimiter='\t');w.writerow(['taxid','reads','low_complexity_reads','bases','dust_bases'])
  for tax,c in taxstats.items():w.writerow([tax,c['reads'],c['low_complexity_reads'],c['bases'],c['dust_bases']])
 result={'source_id':sid,'library_id':row['library_id'],'library_type':row['library_type'],
   'run_id':row['run_id'],'counts':dict(counts),'method':{'dust_window':64,'dust_level':20,
   'flag':'DUST fraction>=0.50 OR Shannon entropy<1.2 bits OR dominant base>=0.80 OR non-DUST bases<100',
   'flags_are_heuristics_not_organism_diagnosis':True}}
 done.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(process,rows))
(OUT/'complexity_summary.json').write_text(json.dumps(results,indent=2)+'\n')
print('All dehosted reads audited',flush=True)
