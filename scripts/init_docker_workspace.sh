#!/usr/bin/env bash
set -euo pipefail

# Initialize the user's ./data directory structure for Docker deployment.
DATA_ROOT="${1:-./data}"

mkdir -p \
  "$DATA_ROOT/input" \
  "$DATA_ROOT/output" \
  "$DATA_ROOT/databases/kraken2" \
  "$DATA_ROOT/databases/viral_refseq_prot" \
  "$DATA_ROOT/references/pathogen_panel"

echo "Initialized Docker data workspace under: $(cd "$DATA_ROOT" && pwd)"
echo "  - $DATA_ROOT/input       (原始樣本：<sample_id>/dna/*.fastq.gz, <sample_id>/cdna/*.fastq.gz)"
echo "  - $DATA_ROOT/output      (分析產出結果與 jobs.db)"
echo "  - $DATA_ROOT/databases   (Kraken2 PlusPF 與 viral_refseq_prot)"
echo "  - $DATA_ROOT/references  (human-masked 與 pathogen_panel)"

