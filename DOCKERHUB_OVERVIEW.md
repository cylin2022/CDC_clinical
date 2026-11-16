# ClinMeta — Offline Clinical mNGS Pathogen Analysis

**ClinMeta** is a **100% offline-capable** Nanopore mNGS pathogen analysis platform with a built-in **Web GUI** and **multi-sample FIFO queue**.

---

## Pull & Run (Quick Start)

### 1. Pull the Image & Create `./data/`

```bash
docker pull <your-dockerhub-user>/clinmeta:latest
mkdir -p data/{input,output,databases,references}
```

### 2. Start the Container

**Option A — Standard Offline Run** (mount existing local databases under `./data/`):

```bash
docker run -d \
  --name clinmeta-offline \
  -p 8080:8080 \
  -v "$PWD/data/databases:/data/databases" \
  -v "$PWD/data/references:/data/references" \
  -v "$PWD/data/input:/data/input:ro" \
  -v "$PWD/data/output:/data/output" \
  <your-dockerhub-user>/clinmeta:latest
```

**Option B — First-Time Auto-Download** (downloads reference bundle from Zenodo + `Kraken2 PlusPF` from AWS S3 on startup):

```bash
docker run -d \
  --name clinmeta-offline \
  -p 8080:8080 \
  -e CDC_ZENODO_BUNDLE_URL="<your-zenodo-bundle-url>" \
  -e CDC_AUTO_DOWNLOAD_K2=true \
  -v "$PWD/data/databases:/data/databases" \
  -v "$PWD/data/references:/data/references" \
  -v "$PWD/data/input:/data/input:ro" \
  -v "$PWD/data/output:/data/output" \
  <your-dockerhub-user>/clinmeta:latest
```

Open **`http://<server-ip>:8080`** in your browser to queue samples and view results.

---

## Directory Layout (`./data/`)

```text
./data/
  ├── input/       # <sample_id>/dna/*.fastq.gz & <sample_id>/cdna/*.fastq.gz
  ├── output/      # Analysis results, BAMs, HTML/Sankey reports, jobs.db
  ├── databases/   # kraken2/k2_pluspf_20261022/ (~117 GB) + viral_refseq_prot/ (~331 MB)
  └── references/  # human-masked* (~20 GB) + pathogen_panel/ (~136 MB)
```

---

## Pipeline & Included Tools

- **Workflow**: `chopper` QC → `hostile` host removal → `Kraken2` → competitive `minimap2` validation → `DIAMOND blastx` unknown-virus screen → **Sankey & HTML reports**
- **Tools**: `Python 3.12` · `chopper 0.15.0` · `hostile 2.0.2` · `minimap2 2.31` · `kraken2 2.17.2` · `samtools 1.25` · `diamond 2.1.9` · `seqkit 2.14` · `pigz 2.8`

---

## Environment Variables & Requirements

| Variable | Default | Purpose |
|---|---|---|
| `CDC_PORT` | `8080` | Web GUI port |
| `CDC_INPUT_DIR` / `CDC_OUTPUT_DIR` | `/data/input` / `/data/output` | Sample I/O paths |
| `CDC_DB_DIR` / `CDC_REF_DIR` | `/data/databases` / `/data/references` | Database & reference paths |
| `CDC_ZENODO_BUNDLE_URL` | `""` | Auto-download ref bundle if missing |
| `CDC_AUTO_DOWNLOAD_K2` | `false` | Auto-download Kraken2 PlusPF from AWS S3 |
| `CDC_BUILD_THREADS` | `16` | Threads for startup index building |

- **Hardware**: **≥ 128 GB RAM** (for Kraken2 PlusPF), **≥ 200 GB SSD/NVMe**, **16+ CPU cores**
