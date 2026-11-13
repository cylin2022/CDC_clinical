#!/usr/bin/env python3
"""Select candidate species reproducibly and retrieve database-matched public genomes."""
import csv,gzip,hashlib,json,re,subprocess,concurrent.futures
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/pathogen_validation'
REF=OUT/'references'
REF.mkdir(parents=True,exist_ok=True)
DB=ROOT/'databases/kraken2/k2_pluspf_20261022'
parents={};rank={}
for line in (DB/'nodes.dmp').open():
 f=line.split('|');parents[int(f[0])]=int(f[1]);rank[int(f[0])]=f[2].strip()
def ancestor(t,target=None):
 seen=set()
 while t not in seen:
  if (target is not None and t==target) or (target is None and rank.get(t)=='species'):return t
  seen.add(t);t=parents.get(t,1)
 return None
rows=list(csv.DictReader((ROOT/'clean_results/kraken2/species_table.tsv').open(),delimiter='\t'))
prefixes=['library.sample_001.SG26-GS01','library.sample_001.ST26-UM01']
selected=[]
for r in rows:
 t=int(r['taxid']);n=max(int(r[p+'.reads']) for p in prefixes);d=max(int(r[p+'.distinct_minimizers']) for p in prefixes)
 if t!=9606 and (n>=500 or d>=500 or (n>=100 and d>=50) or (ancestor(t,10239) and n>=3)):
  r['selection_reason']='reads>=500 OR distinct_minimizers>=500 OR (reads>=100 AND distinct>=50) OR virus_reads>=3'
  selected.append(r)
with (OUT/'candidates.tsv').open('w') as handle:
 w=csv.DictWriter(handle,fieldnames=list(selected[0]),delimiter='\t');w.writeheader();w.writerows(selected)
chosen={int(r['taxid']):r for r in selected}
accessions={}
for line in (DB/'seqid2taxid.map').open():
 seq,t=line.rstrip().split('\t');species=ancestor(int(t))
 if species in chosen:
  accessions[seq.split('|')[-1]]=(species,int(t))
urls={t:{} for t in chosen}
for line in (DB/'library_report.tsv').open():
 fields=line.rstrip('\n').split('\t')
 if len(fields)<3:continue
 accession=fields[1].split(' ')[0]
 if accession not in accessions:continue
 t,strain=accessions[accession];url=fields[2]
 if not url.startswith('https://ftp.ncbi.nlm.nih.gov/genomes/all/') or not url.endswith('_genomic.fna.gz'):continue
 urls[t].setdefault(url,[]).append((accession,strain,fields[1]))
# One database-matched assembly per species; all available Influenza A reference assemblies.
# Add distinct reference assembly from NCBI Datasets for well-supported priority species.
jobs=[];missing=[]
for t,r in chosen.items():
 entries=list(urls[t].items())
 if not entries:
  missing.append({'taxid':t,'name':r['name'],'reason':'No assembly URL in database library report'})
  continue
 if t==2955291:
  picked=entries
 else:
  # Prefer assemblies containing more of the database's candidate sequences.
  picked=[max(entries,key=lambda x:len(x[1]))]
 for url,seqs in picked:
  assembly=url.split('/')[-2]
  jobs.append({'taxid':t,'name':r['name'],'assembly':assembly,'url':url,
               'database_sequences':len(seqs),'accessions':[s[0] for s in seqs],
               'path':str(REF/(str(t)+'__'+assembly+'.fna.gz'))})
def fetch(job):
 p=Path(job['path'])
 if not p.exists():
  temp=p.with_suffix('.download')
  subprocess.run(['curl','-fsSL','--retry','4','--connect-timeout','30',job['url'],'-o',str(temp)],check=True)
  with gzip.open(temp,'rb') as handle:
   if not handle.read(1)==b'>':raise ValueError('Not FASTA: '+job['url'])
  temp.replace(p)
 with p.open('rb') as handle:job['sha256']=hashlib.file_digest(handle,'sha256').hexdigest()
 print('ready',job['taxid'],job['assembly'],flush=True)
 return job
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for job in pool.map(fetch,jobs):results.append(job)
(OUT/'reference_manifest.json').write_text(json.dumps({'selection_count':len(selected),
 'references':results,'missing':missing,'reference_choice':'database-matched assembly with most sequence entries; all Influenza A assemblies'},indent=2)+'\n')
print('Species:',len(selected),'assemblies:',len(results),'missing:',missing,flush=True)
# Stable contig labels prevent collisions between assemblies; metadata preserves original accession.
contigs=[]
with (REF/'microbes.fa').open('w') as out:
 for job in results:
  with gzip.open(job['path'],'rt') as handle:
   for line in handle:
    if line.startswith('>'):
     original=line[1:].strip();acc=original.split()[0]
     label=f"tax{job['taxid']}__{job['assembly']}__{acc}"
     contigs.append({'contig':label,'taxid':job['taxid'],'species':job['name'],
                     'assembly':job['assembly'],'original_header':original})
     out.write('>'+label+'\n')
    else:out.write(line.upper())
with (REF/'contigs.tsv').open('w') as handle:
 w=csv.DictWriter(handle,fieldnames=list(contigs[0]),delimiter='\t');w.writeheader();w.writerows(contigs)
print('Built microbial reference panel',flush=True)
