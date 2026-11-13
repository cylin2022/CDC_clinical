#!/usr/bin/env python3
"""Integration control: host removal, microbial retention, short/low-Q filtering."""
import csv
import gzip
import hashlib
import urllib.request
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import nanopore_pipeline as p
work=ROOT/'validation/controls'
work.mkdir(parents=True,exist_ok=True)
parts=[]
length=0
with gzip.open(ROOT/'references/human-masked.fa.gz','rt') as handle:
    first=handle.readline().strip()
    for line in handle:
        if line.startswith('>'):
            break
        parts.append(line.strip().upper())
        length+=len(parts[-1])
        if length >= 2_000_000:
            break
seq=''.join(parts)
region=re.search('[ACGT]{10000}',seq[1_000_000:])
if region is None:
    raise SystemExit('Cannot find suitable synthetic human control region')
human=region.group()
donor=human.index('GT',1000)
acceptor=human.index('AG',donor+2000)+2
spliced=human[donor-500:donor]+human[acceptor:acceptor+500]
lambda_path=ROOT/'references/lambda-control.fa'
if not lambda_path.exists():
    urllib.request.urlretrieve('https://www.ebi.ac.uk/ena/browser/api/fasta/J02459.1?download=false',lambda_path)
lambda_seq=''.join(line.strip() for line in (ROOT/'references/lambda-control.fa').read_text().splitlines()
                   if not line.startswith('>')).upper()
if len(lambda_seq) != 48502 or hashlib.sha256(lambda_seq.encode()).hexdigest() != '36432a40f602258d19ae7c8152ddbc30390b559f2859c01d7047c77b048c71b3':
    raise SystemExit('Unexpected lambda control reference')
records=[('human_genomic',human[100:1100],'I'),('human_spliced',spliced,'I'),
         ('lambda_retained',lambda_seq[1000:2000],'I'),
         ('low_quality',lambda_seq[3000:4000],'!'),('too_short',lambda_seq[5000:5100],'I')]
rows=[]
for barcode,library,kind in [('barcode17','SG26-GS01','DNA'),('barcode18','ST26-UM01','cDNA')]:
    source=work/f'{kind}.fastq.gz'
    content=''.join(f'@{name}\n{sequence}\n+\n{quality*len(sequence)}\n'
                    for name,sequence,quality in records)
    if source.exists():
        with gzip.open(source,'rt') as handle:
            if handle.read() != content:
                raise SystemExit('Control input changed; use a fresh validation directory')
    else:
        with gzip.open(source,'wt') as out:
            out.write(content)
    rows.append(dict(sample_id='control',source_id=kind,run_id='synthetic',barcode=barcode,
                     library_id=library,library_type=kind,fastq_path=str(source),size_bytes=source.stat().st_size))
manifest=work/'samples.tsv'
p.write_tsv(manifest,rows,p.FIELDS)
subprocess.run(['bash',str(ROOT/'batch_dehost.sh'),'--manifest',str(manifest),
                '--output',str(work/'results'),'--threads','4','--keep-qc'],check=True)
with (work/'results/source_summary.tsv').open() as handle:
    for row in csv.DictReader(handle,delimiter='\t'):
        assert (int(row['raw_reads']),int(row['qc_reads']),int(row['nonhost_reads']),
                int(row['host_removed_reads'])) == (5,3,1,2),row
        with gzip.open(work/'results/sources'/row['source_id']/'clean.fastq.gz','rt') as reads:
            assert reads.readline().strip() == '@lambda_retained'
p.write_json(work/'validation.json',{'passed':True,'human_reference_first_contig':first,
    'lambda_accession':'J02459.1','expected_per_library':{'raw_reads':5,'qc_reads':3,
    'host_removed_reads':2,'nonhost_reads':1},'retained_read':'lambda_retained'})
print('PASS: DNA and cDNA human controls removed, lambda retained, low-Q and short reads filtered')
