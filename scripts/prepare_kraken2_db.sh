#!/usr/bin/env bash
# Verify (size + S3 multipart ETag) and extract the Kraken2 PlusPF database.
# Usage: bash scripts/prepare_kraken2_db.sh [wait_for_pid]
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
D="$ROOT/databases/kraken2"
NAME=k2_pluspf_20261022
TAR="$D/$NAME.tar.gz"
HDR="$D/$NAME.headers.txt"
if [[ -n "${1:-}" ]]; then
  while kill -0 "$1" 2>/dev/null; do sleep 60; done
fi
expected_size=$(awk -F': ' 'tolower($1)=="content-length"{print $2}' "$HDR")
etag=$(awk -F': ' 'tolower($1)=="etag"{gsub(/"/,"",$2); print $2}' "$HDR")
actual_size=$(stat -c %s "$TAR")
echo "size expected=$expected_size actual=$actual_size"
[[ "$expected_size" == "$actual_size" ]] || { echo "ERROR: size mismatch" >&2; exit 1; }
python3 - "$TAR" "$etag" <<'EOF'
import hashlib, sys
path, etag = sys.argv[1], sys.argv[2]
want, parts = etag.split('-'); parts = int(parts)
size = 16 * 1024 * 1024
digests = []
with open(path, 'rb') as fh:
    while chunk := fh.read(size):
        digests.append(hashlib.md5(chunk).digest())
got = hashlib.md5(b''.join(digests)).hexdigest()
print(f'ETag expected={want}-{parts} computed={got}-{len(digests)}')
if got != want or len(digests) != parts:
    sys.exit('ERROR: ETag mismatch')
EOF
mkdir -p "$D/$NAME.partial"
tar -I pigz -xf "$TAR" -C "$D/$NAME.partial"
# Tarball may contain files at top level or inside one directory.
if [[ -f "$D/$NAME.partial/hash.k2d" ]]; then
  mv "$D/$NAME.partial" "$D/$NAME"
else
  sub=$(dirname "$(find "$D/$NAME.partial" -name hash.k2d | head -1)")
  mv "$sub" "$D/$NAME"; rm -rf "$D/$NAME.partial"
fi
ls -la "$D/$NAME"
"$ROOT/.conda/envs/nanopore/bin/kraken2-inspect" --db "$D/$NAME" --skip-counts | head -20 \
  > "$D/$NAME.inspect-head.txt" 2>&1 || true
sha256sum "$TAR" > "$TAR.sha256"
echo "READY: $D/$NAME"

