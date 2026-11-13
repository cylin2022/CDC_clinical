#!/usr/bin/env python3
"""Local Nanopore QC/dehosting with source provenance and resumable outputs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import fcntl
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ['sample_id', 'source_id', 'run_id', 'barcode', 'library_id', 'library_type', 'fastq_path', 'size_bytes']


def write_json(path, obj):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')
    temp.replace(path)


def write_tsv(path, rows, fields):
    temp = path.with_suffix('.tsv.tmp')
    with temp.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def scan(args):
    source = args.input.resolve()
    if not source.is_dir():
        raise ValueError(f'Input directory does not exist: {source}')
    files = sorted(p for p in source.rglob('*') if p.is_file() and
                   p.name.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')))
    if not files:
        raise ValueError(f'No FASTQ files found under {source}')
    rows = []
    for path in files:
        stem = re.sub(r'\.(fastq|fq)(\.gz)?$', '', path.name)
        suffix = hashlib.sha256(str(path.relative_to(source)).encode()).hexdigest()[:8]
        library = {'barcode17': ('SG26-GS01', 'DNA'), 'barcode18': ('ST26-UM01', 'cDNA')}
        if path.parent.name not in library:
            raise ValueError(f'Unknown barcode/library: {path.parent.name}')
        library_id, library_type = library[path.parent.name]
        rows.append(dict(sample_id=args.sample_id, source_id=f'{stem}_{suffix}',
                         run_id=path.name.split('_')[0], barcode=path.parent.name,
                         library_id=library_id, library_type=library_type,
                         fastq_path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT)
                         else str(path), size_bytes=path.stat().st_size))
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    if args.manifest.exists():
        raise ValueError(f'Refusing to overwrite existing manifest: {args.manifest}')
    write_tsv(args.manifest, rows, FIELDS)
    print(f'Manifest: {args.manifest}; {len(rows)} sources, one sample: {args.sample_id}')


def load_manifest(path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != FIELDS:
            raise ValueError(f'Manifest columns must be: {FIELDS}')
        rows = list(reader)
    if not rows:
        raise ValueError('Manifest has no sources')
    seen_ids, seen_paths = set(), set()
    for row in rows:
        for field in ('sample_id', 'source_id', 'library_id', 'library_type'):
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', row[field]):
                raise ValueError(f'Unsafe identifier: {row[field]}')
        source = (ROOT / row['fastq_path']).resolve()
        if row['source_id'] in seen_ids or source in seen_paths:
            raise ValueError('Duplicate source identifier or FASTQ path')
        seen_ids.add(row['source_id'])
        seen_paths.add(source)
        if not source.is_file() or source.stat().st_size != int(row['size_bytes']):
            raise ValueError(f'Missing source or changed size: {source}')
        expected = {'barcode17': ('SG26-GS01', 'DNA'), 'barcode18': ('ST26-UM01', 'cDNA')}
        if expected.get(row['barcode']) != (row['library_id'], row['library_type']):
            raise ValueError('Barcode/library mapping does not match confirmed metadata')
        row['path'] = source
    if len({row['sample_id'] for row in rows}) != 1:
        raise ValueError('This workflow expects exactly one biological sample')
    return rows


def subset(source, target, count):
    opener = gzip.open if source.name.endswith('.gz') else open
    temp = target.with_suffix('.tmp')
    reads = 0
    with opener(source, 'rb') as inp, gzip.open(temp, 'wb') as out:
        for _ in range(count):
            header = inp.readline()
            if not header:
                break
            seq, plus, qual = (inp.readline() for _ in range(3))
            if not header.startswith(b'@') or not plus.startswith(b'+') or not qual or \
                    len(seq.rstrip()) != len(qual.rstrip()):
                raise ValueError(f'Invalid four-line FASTQ record in {source}')
            out.writelines((header, seq, plus, qual))
            reads += 1
    if not reads:
        raise ValueError(f'Empty source: {source}')
    temp.replace(target)


def stats(path, threads):
    result = subprocess.run(['seqkit', 'stats', '-T', '-j', str(threads), str(path)],
                            capture_output=True, text=True, check=True)
    row = next(csv.DictReader(result.stdout.splitlines(), delimiter='\t'))
    return {'reads': int(row['num_seqs'].replace(',', '')),
            'bases': int(row['sum_len'].replace(',', ''))}


def qc(source, target, args, log):
    temp = target.with_suffix('.tmp')
    command = ['chopper', '-q', str(args.quality), '-l', str(args.min_length),
               '--threads', str(args.threads)]
    if args.headcrop or args.tailcrop:
        help_text = subprocess.run(['chopper', '--help'], capture_output=True,
                                   text=True, check=True).stdout
        if '--trim-approach' in help_text:
            command += ['--trim-approach', 'fixed-crop']
        command += ['--headcrop', str(args.headcrop), '--tailcrop', str(args.tailcrop)]
    decompress = ['pigz', '-dc', str(source)] if source.name.endswith('.gz') \
        else ['cat', str(source)]
    with log.open('wb') as err, temp.open('wb') as out:
        first = subprocess.Popen(decompress, stdout=subprocess.PIPE, stderr=err)
        second = subprocess.Popen(command, stdin=first.stdout, stdout=subprocess.PIPE, stderr=err)
        first.stdout.close()
        third = subprocess.Popen(['pigz', '-p', str(min(args.threads, 4)), '-c'],
                                 stdin=second.stdout, stdout=out, stderr=err)
        second.stdout.close()
        codes = [third.wait(), second.wait(), first.wait()]
    if any(codes):
        raise RuntimeError(f'QC failed ({codes}); see {log}')
    temp.replace(target)


def run(args):
    rows = load_manifest(args.manifest)
    for tool in ('chopper', 'hostile', 'pigz', 'seqkit', 'minimap2', 'samtools'):
        if not shutil.which(tool):
            raise ValueError(f'Required executable not found: {tool}')
    index = args.index.resolve()
    if not index.is_file():
        raise ValueError(f'Host reference/index not found: {index}')
    if any(row['library_type'] == 'cDNA' for row in rows) and not args.cdna_index.is_file():
        raise ValueError(f'cDNA splice index not found: {args.cdna_index}')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (output / '.pipeline.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        execute(args, rows, index, output)


def execute(args, rows, index, output):
    summary, cleaned = [], {}
    config = {key: getattr(args, key) for key in
              ('threads', 'quality', 'min_length', 'headcrop', 'tailcrop', 'limit_reads')}
    config.update(index=str(index), index_size=index.stat().st_size,
                  index_mtime_ns=index.stat().st_mtime_ns,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if any(row['library_type'] == 'cDNA' for row in rows):
        cdna_index = args.cdna_index.resolve()
        config.update(cdna_index=str(cdna_index), cdna_index_size=cdna_index.stat().st_size,
                      cdna_index_mtime_ns=cdna_index.stat().st_mtime_ns)
    env_lock = ROOT / 'environment-linux-64.lock.txt'
    if env_lock.is_file():
        config['environment_lock_sha256'] = hashlib.sha256(env_lock.read_bytes()).hexdigest()
    config['tools'] = {name: str(Path(shutil.which(name)).resolve()) for name in
                       ('chopper', 'hostile', 'seqkit', 'minimap2', 'samtools', 'pigz')}
    def process(row):
        sid, source = row['source_id'], row['path']
        selected_index = args.cdna_index.resolve() if row['library_type'] == 'cDNA' else index
        preset = 'splice' if row['library_type'] == 'cDNA' else 'map-ont'
        work = output / 'sources' / sid
        work.mkdir(parents=True, exist_ok=True)
        fingerprint = dict(config, source=str(source), size=source.stat().st_size,
                           mtime_ns=source.stat().st_mtime_ns, sample_id=row['sample_id'],
                           library_id=row['library_id'], library_type=row['library_type'])
        done = work / 'complete.json'
        clean = work / 'clean.fastq.gz'
        if done.exists():
            previous = json.loads(done.read_text())
            if previous['fingerprint'] != fingerprint:
                raise ValueError(f'Configuration/source changed for {sid}; use a new output directory')
            if not clean.is_file() or clean.stat().st_size != previous['clean_bytes']:
                raise ValueError(f'Completed output missing or changed for {sid}')
            item = previous['summary']
            print(f'[resume] {sid}', flush=True)
        else:
            print(f'[processing] {sid}', flush=True)
            selected = source
            if args.limit_reads:
                selected = work / 'subset.fastq.gz'
                subset(source, selected, args.limit_reads)
            raw_stats = stats(selected, args.threads)
            qc_path = work / 'qc.fastq.gz'
            qc(selected, qc_path, args, work / 'chopper.log')
            qc_stats = stats(qc_path, args.threads)
            if qc_stats['reads']:
                host_out = work / 'hostile'
                host_out.mkdir(exist_ok=True)
                with (work / 'hostile.json').open('w') as stdout, \
                     (work / 'hostile.log').open('w') as stderr:
                    subprocess.run(['hostile', 'clean', '--fastq1', str(qc_path),
                                    '--aligner', 'minimap2', '--index', str(selected_index),
                                    '--aligner-args=-x ' + preset, '--airplane',
                                    '--threads', str(args.threads), '--output', str(host_out),
                                    '--force'], stdout=stdout, stderr=stderr, check=True)
                candidates = list(host_out.glob('*.fastq.gz'))
                if len(candidates) != 1:
                    raise RuntimeError(f'Expected one Hostile output in {host_out}')
                candidates[0].replace(clean)
            else:
                with gzip.open(clean, 'wb'):
                    pass
                write_json(work / 'hostile.json', {'skipped': 'no reads passed QC'})
            clean_stats = stats(clean, args.threads)
            if not 0 <= clean_stats['reads'] <= qc_stats['reads'] <= raw_stats['reads']:
                raise RuntimeError(f'Inconsistent read counts: {sid}')
            item = dict(sample_id=row['sample_id'], source_id=sid, run_id=row['run_id'],
                        barcode=row['barcode'], library_id=row['library_id'],
                        library_type=row['library_type'], raw_reads=raw_stats['reads'],
                        raw_bases=raw_stats['bases'], qc_reads=qc_stats['reads'],
                        qc_bases=qc_stats['bases'], nonhost_reads=clean_stats['reads'],
                        nonhost_bases=clean_stats['bases'],
                        host_removed_reads=qc_stats['reads'] - clean_stats['reads'])
            write_json(done, dict(fingerprint=fingerprint, summary=item,
                                 clean_bytes=clean.stat().st_size))
            if not args.keep_qc:
                qc_path.unlink()
            print(f'[done] {sid}', flush=True)
        return item, clean

    # Largest sources first so the long DNA runs start immediately.
    order = sorted(rows, key=lambda r: r['path'].stat().st_size, reverse=True)
    with ThreadPoolExecutor(max_workers=max(1, min(args.jobs, len(rows)))) as pool:
        futures = {id(row): pool.submit(process, row) for row in order}
        results = {key: future.result() for key, future in futures.items()}
    for row in rows:  # manifest order for deterministic summaries/merges
        item, clean = results[id(row)]
        summary.append(item)
        cleaned.setdefault(row['library_id'], []).append(clean)
    write_tsv(output / 'source_summary.tsv', summary, list(summary[0]))
    sample_id = rows[0]['sample_id']
    merged_libraries = {}
    for library_id, paths in cleaned.items():
        merged = output / f'{sample_id}.{library_id}.nonhost.fastq.gz'
        temp = merged.with_suffix('.tmp')
        with temp.open('wb') as out:
            for path in paths:
                with path.open('rb') as inp:
                    shutil.copyfileobj(inp, out)
        merged_stats = stats(temp, args.threads)
        selected = [x for x in summary if x['library_id'] == library_id]
        if merged_stats['reads'] != sum(x['nonhost_reads'] for x in selected) or \
                merged_stats['bases'] != sum(x['nonhost_bases'] for x in selected):
            raise RuntimeError('Merged output does not match source totals')
        temp.replace(merged)
        merged_libraries[library_id] = dict(fastq=str(merged), **merged_stats)
    totals = {key: sum(x[key] for x in summary) for key in summary[0]
              if key.endswith(('_reads', '_bases'))}
    write_json(output / 'sample_summary.json', dict(sample_id=sample_id,
               mode='prefix-smoke-test' if args.limit_reads else 'full',
               source_count=len(rows), totals=totals, merged_libraries=merged_libraries, config=config))
    print(json.dumps(totals, indent=2), flush=True)
    print(f'Complete: {output}; libraries kept separate', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    scan_parser = commands.add_parser('scan')
    scan_parser.add_argument('--input', type=Path, default=ROOT / 'fastq_pass')
    scan_parser.add_argument('--sample-id', default='sample_001')
    scan_parser.add_argument('--manifest', type=Path, default=ROOT / 'samples.tsv')
    scan_parser.set_defaults(func=scan)
    runner = commands.add_parser('run')
    runner.add_argument('--manifest', type=Path, default=ROOT / 'samples.tsv')
    runner.add_argument('--output', type=Path, default=ROOT / 'clean_results')
    runner.add_argument('--index', type=Path, default=ROOT / 'references/human-masked.mmi')
    runner.add_argument('--cdna-index', type=Path, default=ROOT / 'references/human-masked.splice.mmi')
    runner.add_argument('--threads', type=int, default=16, help='Threads per source job')
    runner.add_argument('--jobs', type=int, default=1, help='Sources processed in parallel')
    runner.add_argument('--quality', type=float, default=10)
    runner.add_argument('--min-length', type=int, default=200)
    runner.add_argument('--headcrop', type=int, default=15)
    runner.add_argument('--tailcrop', type=int, default=15)
    runner.add_argument('--limit-reads', type=int, default=0,
                        help='Take the first N reads per source; 0 processes full files')
    runner.add_argument('--keep-qc', action='store_true')
    runner.set_defaults(func=run)
    args = parser.parse_args()
    if args.command == 'run' and (args.threads < 1 or args.jobs < 1 or args.min_length < 1 or
            args.limit_reads < 0 or args.headcrop < 0 or args.tailcrop < 0 or args.quality < 0):
        parser.error('Invalid thread count or QC/subset parameter')
    try:
        args.func(args)
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
