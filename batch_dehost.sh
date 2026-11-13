#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$ROOT/.conda/envs/nanopore/bin:$PATH"
export XDG_CACHE_HOME="$ROOT/.cache"
export HOSTILE_CACHE_DIR="$ROOT/.cache/hostile"
exec python "$ROOT/scripts/nanopore_pipeline.py" run "$@"
