<#!/usr/bin/env bash
# Kraken2 classification of dehosted reads: per sequencing source and per merged library.
# Usage: bash scripts/run_kraken2.sh [clean_results] [databases/kraken2/k2_pluspf_20261022] [threads]
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export PATH="$ROOT/.conda/envs/nanopore/bin:$PATH"
RESULTS="$(realpath "${1:-$ROOT/clean_results}")"
DB="$(realpath "${2:-$ROOT/databases/kraken2/k2_pluspf_20261022}")"
THREADS="${3:-64}"
OUT="$RESULTS/kraken2"
for f in hash.k2d opts.k2d taxo.k2d; do
  [[ -s "$DB/$f" ]] || { echo "ERROR: $DB/$f missing" >&2; exit 1; }
done
[[ -s "$RESULTS/sample_summary.json" ]] || { echo "ERROR: dehosting not complete in $RESULTS" >&2; exit 1; }
mkdir -p "$OUT"
kraken2 --version | head -1 > "$OUT/kraken2_version.txt"
echo "$DB" > "$OUT/database_path.txt"
# Warm the page cache once so memory-mapped lookups are served from RAM for every run.
cat "$DB"/*.k2d > /dev/null

classify() {  # name fastq
  local name="$1" fq="$2"
  if [[ -s "$OUT/$name.report.txt" && -s "$OUT/$name.done" ]]; then
    echo "[resume] $name"; return
  fi
  echo "[kraken2] $name"
  kraken2 --db "$DB" --threads "$THREADS" --memory-mapping \
    --report "$OUT/$name.report.tmp" --report-minimizer-data \
    --output >(pigz -p 8 > "$OUT/$name.reads.kraken.gz") \
    --use-names "$fq" 2> "$OUT/$name.kraken2.log"
  mv "$OUT/$name.report.tmp" "$OUT/$name.report.txt"
  date -Is > "$OUT/$name.done"
}

for src in "$RESULTS"/sources/*/; do
  classify "source.$(basename "$src")" "$src/clean.fastq.gz"
done
for fq in "$RESULTS"/*.nonhost.fastq.gz; do
  name="$(basename "$fq" .nonhost.fastq.gz)"
  classify "library.$name" "$fq"
done
python "$ROOT/scripts/summarize_kraken2.py" "$OUT"
python "$ROOT/scripts/kraken2_report_html.py" "$OUT" "$OUT/report.html"
echo "Complete: $OUT (HTML report: $OUT/report.html)"
