#!/usr/bin/env python3
"""Fetch and checksum the official masked human reference and build both presets."""
import hashlib
import json
from pathlib import Path
import subprocess
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
NAME = 'human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401'
URL = 'https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o'
refdir = ROOT / 'references'
refdir.mkdir(exist_ok=True)
(ROOT / 'logs').mkdir(exist_ok=True)
manifest = refdir / 'hostile-index-manifest.json'
if not manifest.exists():
    urllib.request.urlretrieve(URL + '/manifest.json', manifest)
expected = json.loads(manifest.read_text())[NAME]['assets'][NAME + '.fa.gz']['sha256']
reference = refdir / 'human-masked.fa.gz'
if not reference.exists():
    temp = reference.with_suffix('.download')
    subprocess.run(['curl', '-fL', '--retry', '5', '-C', '-', URL + '/' + NAME + '.fa.gz',
                    '-o', str(temp)], check=True)
    temp.replace(reference)
with reference.open('rb') as handle:
    actual = hashlib.file_digest(handle, 'sha256').hexdigest()
if actual != expected:
    raise SystemExit('Human reference SHA256 mismatch; remove corrupt download and retry')
for preset, filename in [('map-ont', 'human-masked.mmi'), ('splice', 'human-masked.splice.mmi')]:
    index = refdir / filename
    if not index.exists():
        temp = index.with_suffix('.mmi.tmp')
        with (ROOT / 'logs' / f'index-{preset}.log').open('w') as log:
            subprocess.run(['minimap2', '-t', '16', '-x', preset, '-d', str(temp),
                            str(reference)], stdout=log, stderr=log, check=True)
        temp.replace(index)
(refdir / 'reference.json').write_text(json.dumps(dict(index_name=NAME, sha256=actual,
    source_url=URL + '/' + NAME + '.fa.gz', minimap2_presets=['map-ont', 'splice']), indent=2)+'\n')
print('Human reference verified; DNA and cDNA indexes ready')
