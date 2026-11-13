#!/usr/bin/env bash
set -euo pipefail

REF_DIR="${CDC_REF_DIR:-/data/references}"
DB_DIR="${CDC_DB_DIR:-/data/databases}"
ZENODO_URL="${CDC_ZENODO_BUNDLE_URL:-}"
ZENODO_SHA256="${CDC_ZENODO_BUNDLE_SHA256:-}"
AUTO_DL_K2="${CDC_AUTO_DOWNLOAD_K2:-false}"
THREADS="${CDC_BUILD_THREADS:-16}"

INPUT_DIR="${CDC_INPUT_DIR:-/data/input}"
OUTPUT_DIR="${CDC_OUTPUT_DIR:-/data/output}"

mkdir -p "$INPUT_DIR" "$OUTPUT_DIR" "$REF_DIR" "$DB_DIR/kraken2" "$DB_DIR/viral_refseq_prot"

echo "[bootstrap] Checking lightweight reference & viral protein bundle..."
if [[ ! -f "$REF_DIR/human-masked.fa.gz" || ! -f "$DB_DIR/viral_refseq_prot/viral.1.protein.faa.gz" ]]; then
  if [[ -n "$ZENODO_URL" ]]; then
    TMP_TAR="/tmp/zenodo_bundle.tar.gz"
    echo "[bootstrap] Downloading Zenodo bundle from $ZENODO_URL"
    curl -L -o "$TMP_TAR" "$ZENODO_URL"
    if [[ -n "$ZENODO_SHA256" ]]; then
      echo "$ZENODO_SHA256  $TMP_TAR" | sha256sum -c -
    fi
    tar -xzf "$TMP_TAR" -C /tmp/zenodo_extract
    cp -rn /tmp/zenodo_extract/references/* "$REF_DIR/" || true
    cp -rn /tmp/zenodo_extract/databases/* "$DB_DIR/" || true
    rm -rf "$TMP_TAR" /tmp/zenodo_extract
  else
    echo "[bootstrap] WARNING: Reference bundle missing and CDC_ZENODO_BUNDLE_URL not set."
  fi
fi

# Build minimap2 indices locally if human-masked.fa.gz is present (saves ~19 GB transfer)
if [[ -f "$REF_DIR/human-masked.fa.gz" && ! -f "$REF_DIR/human-masked.mmi" ]]; then
  echo "[bootstrap] Building human-masked.mmi (DNA map-ont)..."
  minimap2 -x map-ont -t "$THREADS" -d "$REF_DIR/human-masked.mmi" "$REF_DIR/human-masked.fa.gz"
fi
if [[ -f "$REF_DIR/human-masked.fa.gz" && ! -f "$REF_DIR/human-masked.splice.mmi" ]]; then
  echo "[bootstrap] Building human-masked.splice.mmi (cDNA splice)..."
  minimap2 -x splice -t "$THREADS" -d "$REF_DIR/human-masked.splice.mmi" "$REF_DIR/human-masked.fa.gz"
fi

# Build DIAMOND viral protein DB locally if .faa.gz is present
if [[ -f "$DB_DIR/viral_refseq_prot/viral.1.protein.faa.gz" && ! -f "$DB_DIR/viral_refseq_prot/viral_refseq_prot.dmnd" ]]; then
  echo "[bootstrap] Building viral_refseq_prot.dmnd..."
  diamond makedb --in "$DB_DIR/viral_refseq_prot/viral.1.protein.faa.gz" \
    -d "$DB_DIR/viral_refseq_prot/viral_refseq_prot" --threads "$THREADS"
fi

# Kraken2 PlusPF (117 GB) from AWS S3 (not Zenodo)
K2_HASH="$DB_DIR/kraken2/k2_pluspf_20261022/hash.k2d"
if [[ ! -f "$K2_HASH" ]]; then
  if [[ "$AUTO_DL_K2" == "true" ]]; then
    echo "[bootstrap] Downloading Kraken2 PlusPF (2026-10-22) from AWS S3..."
    bash /app/scripts/prepare_kraken2_db.sh "$DB_DIR/kraken2/k2_pluspf_20261022"
  else
    echo "[bootstrap] NOTE: Kraken2 PlusPF missing at $K2_HASH (set CDC_AUTO_DOWNLOAD_K2=true to fetch from AWS S3)."
  fi
else
  echo "[bootstrap] Kraken2 PlusPF ready."
fi

exec "$@"

