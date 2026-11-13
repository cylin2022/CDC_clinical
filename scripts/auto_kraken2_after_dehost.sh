#!/usr/bin/env bash
# Wait for DB download (pid $1) and dehosting (pid $2), then verify DB and classify.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
bash "$ROOT/scripts/prepare_kraken2_db.sh" "$1"
while kill -0 "$2" 2>/dev/null; do sleep 120; done
grep -q '"mode": "full"' "$ROOT/clean_results/sample_summary.json"
bash "$ROOT/scripts/run_kraken2.sh" "$ROOT/clean_results" "$ROOT/databases/kraken2/k2_pluspf_20261022" 64
