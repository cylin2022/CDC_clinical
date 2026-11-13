#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
CONDA_EXE="${CONDA_EXE:-/home/cylin/miniconda3/bin/conda}"
export CONDA_CHANNEL_PRIORITY=strict
export CONDA_PKGS_DIRS="$ROOT/.conda/pkgs"
mkdir -p logs
if [[ ! -x .conda/envs/nanopore/bin/hostile ]]; then
    "$CONDA_EXE" env create --solver libmamba --prefix "$ROOT/.conda/envs/nanopore" --file environment.yml
fi
export PATH="$ROOT/.conda/envs/nanopore/bin:$PATH"
export HOSTILE_CACHE_DIR="$ROOT/.cache/hostile"
"$CONDA_EXE" list --prefix "$ROOT/.conda/envs/nanopore" --explicit > environment-linux-64.lock.txt
python scripts/prepare_reference.py
