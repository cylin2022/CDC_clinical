#!/usr/bin/env bash
set -euo pipefail

# Package the lightweight (~1.2 GB) reference + viral protein bundle for Zenodo upload.
# Excludes .mmi/.dmnd (built locally on startup) and Kraken2 PlusPF (fetched from AWS S3).

OUT_TAR="${1:-zenodo_cdc_nanopore_refs_v1.tar.gz}"
STAGE_DIR="/tmp/cdc_zenodo_stage"

rm -rf "$STAGE_DIR"
mkdir -p "$STAGE_DIR/references/pathogen_panel" "$STAGE_DIR/databases/viral_refseq_prot"

cp references/human-masked.fa.gz references/hostile-index-manifest.json references/reference.json "$STAGE_DIR/references/"
cp references/pathogen_panel/*.fna.gz references/pathogen_panel/manifest.json "$STAGE_DIR/references/pathogen_panel/"
cp databases/viral_refseq_prot/viral.1.protein.faa.gz "$STAGE_DIR/databases/viral_refseq_prot/"

tar -C "$STAGE_DIR" -czf "$OUT_TAR" references databases
sha256sum "$OUT_TAR" > "${OUT_TAR}.sha256"
rm -rf "$STAGE_DIR"

ls -lh "$OUT_TAR" "${OUT_TAR}.sha256"
cat "${OUT_TAR}.sha256"

