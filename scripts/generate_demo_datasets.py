#!/usr/bin/env python3
"""Generate synthetic Nanopore-like demo FASTQ.GZ datasets with known pathogen spikes.

Outputs into input_samples/<demo_id>/{dna,cdna}/ for immediate Web GUI scanning
and offline pipeline testing, plus input_samples/demo_ground_truth.json.
"""

import gzip
import json
import os
import random

SEED = 42
PANEL_DIR = "references/pathogen_panel"
HUMAN_FA = "references/human-masked.fa.gz"
OUT_ROOT = "input_samples"


def load_fasta_seqs(path, max_contigs=None, min_len=400):
    seqs = []
    cur_name = None
    buf = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_name is not None:
                    s = "".join(buf).upper()
                    if len(s) >= min_len and s.count("N") / len(s) < 0.2:
                        seqs.append((cur_name, s))
                        if max_contigs and len(seqs) >= max_contigs:
                            return seqs
                cur_name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur_name is not None:
            s = "".join(buf).upper()
            if len(s) >= min_len and s.count("N") / len(s) < 0.2:
                seqs.append((cur_name, s))
    return seqs


def mutate_ont(seq, rng, err_rate=0.025):
    """Introduce ~2.5% Nanopore-like substitution/indel noise."""
    bases = ("A", "C", "G", "T")
    out = []
    for ch in seq:
        if ch not in bases:
            out.append(rng.choice(bases))
            continue
        r = rng.random()
        if r < err_rate * 0.6:
            out.append(rng.choice([b for b in bases if b != ch]))
        elif r < err_rate * 0.8:
            continue  # deletion
        elif r < err_rate:
            out.append(ch)
            out.append(rng.choice(bases))  # insertion
        else:
            out.append(ch)
    return "".join(out)


def sample_reads(seqs, n_reads, min_l, max_l, prefix, rng):
    reads = []
    for i in range(n_reads):
        cname, cseq = rng.choice(seqs)
        if len(cseq) <= min_l:
            frag = cseq
        else:
            rlen = rng.randint(min_l, min(max_l, len(cseq)))
            start = rng.randint(0, len(cseq) - rlen)
            frag = cseq[start : start + rlen]
        if rng.random() < 0.5:
            comp = str.maketrans("ACGT", "TGCA")
            frag = frag.translate(comp)[::-1]
        noisy = mutate_ont(frag, rng)
        qual = "".join(chr(rng.randint(47, 56)) for _ in range(len(noisy)))  # Q14-Q23
        reads.append((f"{prefix}_{i:05d} ref={cname}", noisy, qual))
    return reads


def make_tso_concatemer_reads(n_reads, rng):
    motif = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
    reads = []
    for i in range(n_reads):
        seq = (motif * rng.randint(10, 22))[: rng.randint(320, 650)]
        noisy = mutate_ont(seq, rng, err_rate=0.02)
        qual = "".join(chr(rng.randint(47, 56)) for _ in range(len(noisy)))
        reads.append((f"tso_artifact_{i:05d}", noisy, qual))
    return reads


def write_fastq_gz(path, reads):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        for rid, seq, qual in reads:
            f.write(f"@{rid}\n{seq}\n+\n{qual}\n")


def main():
    rng = random.Random(SEED)
    print("Loading human-masked contigs (first 15)...")
    human_seqs = load_fasta_seqs(HUMAN_FA, max_contigs=15, min_len=5000)

    print("Loading pathogen references from panel...")
    h1n1_seqs = load_fasta_seqs(os.path.join(PANEL_DIR, "influenza_a_h1n1_ca07.fna.gz"), min_len=500)
    kp_seqs = load_fasta_seqs(os.path.join(PANEL_DIR, "klebsiella_pneumoniae_hs11286.fna.gz"), min_len=2000)
    candida_seqs = load_fasta_seqs(os.path.join(PANEL_DIR, "candida_albicans_sc5314.fna.gz"), min_len=5000)
    pa_seqs = load_fasta_seqs(os.path.join(PANEL_DIR, "287__GCF_051123625.1_ASM5112362v1.fna.gz"), min_len=2000)

    demos = []

    # 1. demo_flu_h1n1
    dna_1 = sample_reads(human_seqs, 500, 800, 2500, "human_dna", rng)
    cdna_1 = sample_reads(human_seqs, 400, 400, 1400, "human_cdna", rng) + sample_reads(
        h1n1_seqs, 300, 500, 1600, "flu_h1n1", rng
    )
    rng.shuffle(cdna_1)
    write_fastq_gz(os.path.join(OUT_ROOT, "demo_flu_h1n1/dna/demo_flu_h1n1.dna.fastq.gz"), dna_1)
    write_fastq_gz(os.path.join(OUT_ROOT, "demo_flu_h1n1/cdna/demo_flu_h1n1.cdna.fastq.gz"), cdna_1)
    demos.append({
        "sample_id": "demo_flu_h1n1",
        "expected_primary": "Influenza A virus (A/California/07/2009(H1N1))",
        "dna_reads": len(dna_1),
        "cdna_reads": len(cdna_1),
        "spike_details": {"FLU_H1N1_cdna": 300, "human_dna": 500, "human_cdna": 400},
    })

    # 2. demo_kp_bacteremia
    dna_2 = sample_reads(human_seqs, 500, 800, 2500, "human_dna", rng) + sample_reads(
        kp_seqs, 300, 1000, 3500, "kp_dna", rng
    )
    cdna_2 = sample_reads(human_seqs, 400, 400, 1400, "human_cdna", rng) + sample_reads(
        kp_seqs, 150, 400, 1400, "kp_cdna", rng
    )
    rng.shuffle(dna_2)
    rng.shuffle(cdna_2)
    write_fastq_gz(os.path.join(OUT_ROOT, "demo_kp_bacteremia/dna/demo_kp_bacteremia.dna.fastq.gz"), dna_2)
    write_fastq_gz(os.path.join(OUT_ROOT, "demo_kp_bacteremia/cdna/demo_kp_bacteremia.cdna.fastq.gz"), cdna_2)
    demos.append({
        "sample_id": "demo_kp_bacteremia",
        "expected_primary": "Klebsiella pneumoniae",
        "dna_reads": len(dna_2),
        "cdna_reads": len(cdna_2),
        "spike_details": {"Klebsiella_pneumoniae_dna": 300, "Klebsiella_pneumoniae_cdna": 150},
    })

    # 3. demo_candida_coinfection
    dna_3 = (
        sample_reads(human_seqs, 400, 800, 2500, "human_dna", rng)
        + sample_reads(candida_seqs, 200, 1000, 3000, "candida_dna", rng)
        + sample_reads(pa_seqs, 150, 1000, 3000, "paeruginosa_dna", rng)
    )
    cdna_3 = (
        sample_reads(human_seqs, 400, 400, 1400, "human_cdna", rng)
        + sample_reads(candida_seqs, 100, 400, 1200, "candida_cdna", rng)
        + sample_reads(pa_seqs, 100, 400, 1200, "paeruginosa_cdna", rng)
        + make_tso_concatemer_reads(50, rng)
    )
    rng.shuffle(dna_3)
    rng.shuffle(cdna_3)
    write_fastq_gz(os.path.join(OUT_ROOT, "demo_candida_coinfection/dna/demo_candida_coinfection.dna.fastq.gz"), dna_3)
    write_fastq_gz(os.path.join(OUT_ROOT, "demo_candida_coinfection/cdna/demo_candida_coinfection.cdna.fastq.gz"), cdna_3)
    demos.append({
        "sample_id": "demo_candida_coinfection",
        "expected_primary": "Candida albicans + Pseudomonas aeruginosa (+ TSO primer artifacts)",
        "dna_reads": len(dna_3),
        "cdna_reads": len(cdna_3),
        "spike_details": {
            "Candida_albicans_dna": 200,
            "Candida_albicans_cdna": 100,
            "Pseudomonas_aeruginosa_dna": 150,
            "Pseudomonas_aeruginosa_cdna": 100,
            "TSO_primer_artifacts_cdna": 50,
        },
    })

    gt_path = os.path.join(OUT_ROOT, "demo_ground_truth.json")
    with open(gt_path, "w") as f:
        json.dump(demos, f, indent=2, ensure_ascii=False)
    print(f"Generated 3 demo datasets in {OUT_ROOT}/; ground truth -> {gt_path}")


if __name__ == "__main__":
    main()

