#!/usr/bin/env python3
"""Screen for divergent/unknown viruses via:
1) Kraken2 viral Family/Genus rank-level expansion checks (cases where species are fragmented).
2) Primer-filtered DIAMOND blastx against NCBI Viral RefSeq proteins on unmapped/unclassified reads.
"""

import argparse
import collections
import csv
import glob
import gzip
import json
import os
import subprocess


TSO_MOTIFS = (
    "GTATCAACGCAGAGT",
    "ACTCTGCGTTGATAC",
    "AGTGGTATCAAC",
    "GTTGATACCACT",
    "AAGCAGTGGTAT",
    "ATACCACTGCTT",
)


def parse_viral_ranks(report_path):
    """Return dict taxid -> {rank, name, clade_reads, distinct_minimizers} under Viruses (10239)."""
    res = {}
    in_viruses = False
    with open(report_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            clade_reads = int(parts[1])
            dmin = int(parts[4])
            rank = parts[5]
            taxid = parts[6]
            name_raw = parts[7]
            indent = len(name_raw) - len(name_raw.lstrip(" "))
            if taxid == "10239":
                in_viruses = True
                continue
            if in_viruses and indent <= 2:
                break
            if in_viruses and rank in ("F", "G"):
                res[taxid] = {
                    "rank": rank,
                    "name": name_raw.strip(),
                    "clade_reads": clade_reads,
                    "distinct_minimizers": dmin,
                }
    return res


def scan_viral_ranks(k2_dir, out_tsv, min_reads=20, min_dmin=100):
    """Flag viral Families/Genera with >= min_reads and >= min_dmin across libraries/replicates."""
    lib_dna = parse_viral_ranks(os.path.join(k2_dir, "library.sample_001.SG26-GS01.report.txt"))
    lib_cdna = parse_viral_ranks(os.path.join(k2_dir, "library.sample_001.ST26-UM01.report.txt"))

    dna_sources = [parse_viral_ranks(p) for p in sorted(glob.glob(os.path.join(k2_dir, "source.*barcode17*.report.txt")))]
    cdna_sources = [parse_viral_ranks(p) for p in sorted(glob.glob(os.path.join(k2_dir, "source.*barcode18*.report.txt")))]

    all_taxids = set(lib_dna) | set(lib_cdna)
    rows = []
    for tid in all_taxids:
        d = lib_dna.get(tid, {"rank": "", "name": "", "clade_reads": 0, "distinct_minimizers": 0})
        c = lib_cdna.get(tid, {"rank": "", "name": "", "clade_reads": 0, "distinct_minimizers": 0})
        rank = c["rank"] or d["rank"]
        name = c["name"] or d["name"]

        dna_rep = all(s.get(tid, {}).get("clade_reads", 0) >= 3 for s in dna_sources) if dna_sources else False
        cdna_rep = all(s.get(tid, {}).get("clade_reads", 0) >= 3 for s in cdna_sources) if cdna_sources else False

        triggered = (
            (c["clade_reads"] >= min_reads and c["distinct_minimizers"] >= min_dmin)
            or (d["clade_reads"] >= min_reads and d["distinct_minimizers"] >= min_dmin)
        )
        if triggered:
            rows.append({
                "taxid": tid,
                "rank": rank,
                "name": name,
                "cDNA_reads": c["clade_reads"],
                "cDNA_distinct_minimizers": c["distinct_minimizers"],
                "cDNA_both_runs_ge3": "yes" if cdna_rep else "no",
                "DNA_reads": d["clade_reads"],
                "DNA_distinct_minimizers": d["distinct_minimizers"],
                "DNA_both_runs_ge3": "yes" if dna_rep else "no",
                "expansion_recommended": "yes",
            })

    rows.sort(key=lambda r: max(r["cDNA_distinct_minimizers"], r["DNA_distinct_minimizers"]), reverse=True)
    with open(out_tsv, "w", newline="") as f:
        cols = [
            "taxid", "rank", "name",
            "cDNA_reads", "cDNA_distinct_minimizers", "cDNA_both_runs_ge3",
            "DNA_reads", "DNA_distinct_minimizers", "DNA_both_runs_ge3",
            "expansion_recommended",
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


def extract_dna_unclassified(k2_reads_gz, fastq_gz):
    u_ids = set()
    with gzip.open(k2_reads_gz, "rt") as f:
        for line in f:
            if line.startswith("U\t"):
                u_ids.add(line.split("\t", 2)[1])

    seqs = []
    with gzip.open(fastq_gz, "rt") as f:
        while True:
            hdr = f.readline()
            if not hdr:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            rid = hdr[1:].strip().split()[0]
            if rid in u_ids:
                seqs.append((rid, seq))
    return seqs


def read_fasta(path):
    cur = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    yield cur, "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur is not None:
            yield cur, "".join(buf)


def has_tso_motif(seq):
    u = seq.upper()
    return any(m in u for m in TSO_MOTIFS)


def run_diamond_screen(query_fa, dmnd_db, out_tsv, threads):
    raw_hits = out_tsv + ".raw.tsv"
    cmd = [
        "diamond", "blastx",
        "--db", dmnd_db,
        "--query", query_fa,
        "--sensitive",
        "--evalue", "1e-5",
        "--max-target-seqs", "3",
        "--threads", str(threads),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle",
        "--out", raw_hits,
    ]
    subprocess.run(cmd, check=True)

    # Keep best hit per query and filter for meaningful viral protein alignments (aa length >= 50, bitscore >= 50)
    best_by_q = {}
    all_hit_count = 0
    with open(raw_hits) as f:
        for line in f:
            all_hit_count += 1
            cols = line.rstrip("\n").split("\t", 6)
            qid, sid, pident, aln_len, evalue, bitscore, stitle = (
                cols[0], cols[1], float(cols[2]), int(cols[3]), float(cols[4]), float(cols[5]), cols[6]
            )
            if qid not in best_by_q or bitscore > best_by_q[qid]["bitscore"]:
                best_by_q[qid] = {
                    "qseqid": qid,
                    "sseqid": sid,
                    "pident": pident,
                    "aa_length": aln_len,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "stitle": stitle,
                }

    strong_hits = [h for h in best_by_q.values() if h["aa_length"] >= 50 and h["bitscore"] >= 50.0]
    strong_hits.sort(key=lambda h: h["bitscore"], reverse=True)

    with open(out_tsv, "w", newline="") as f:
        cols = ["qseqid", "sseqid", "pident", "aa_length", "evalue", "bitscore", "stitle"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(strong_hits)

    os.remove(raw_hits)
    return len(best_by_q), strong_hits


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--dmnd-db", default="databases/viral_refseq_prot/viral_refseq_prot.dmnd")
    p.add_argument("--outdir", default="analysis/pathogen_validation")
    args = p.parse_args()

    # 1. Viral Family/Genus rank expansion check
    rank_tsv = os.path.join(args.outdir, "viral_rank_expansion.tsv")
    rank_rows = scan_viral_ranks("clean_results/kraken2", rank_tsv)
    print(f"[rank] Flagged {len(rank_rows)} viral Families/Genera -> {rank_tsv}")

    # 2. Prepare primer-filtered unknown FASTA (cDNA truly unmapped K2_U non-TSO + DNA K2_U)
    cdna_kept = []
    cdna_tso = 0
    for rid, seq in read_fasta("analysis/unclassified_check/truly_unmapped.fasta"):
        if has_tso_motif(seq):
            cdna_tso += 1
        else:
            cdna_kept.append((f"cDNA|{rid}", seq))

    dna_u = extract_dna_unclassified(
        "clean_results/kraken2/library.sample_001.SG26-GS01.reads.kraken.gz",
        "clean_results/sample_001.SG26-GS01.nonhost.fastq.gz",
    )
    dna_kept = [(f"DNA|{rid}", seq) for rid, seq in dna_u]

    query_fa = "analysis/unclassified_check/unknown_for_diamond.fasta"
    with open(query_fa, "w") as out:
        for rid, seq in cdna_kept + dna_kept:
            out.write(f">{rid}\n{seq}\n")

    print(f"[query] Wrote {len(cdna_kept)} cDNA + {len(dna_kept)} DNA = {len(cdna_kept)+len(dna_kept)} reads -> {query_fa}")

    # 3. DIAMOND blastx against Viral RefSeq proteins
    hits_tsv = os.path.join(args.outdir, "unknown_viral_screen.tsv")
    raw_q_hits, strong_hits = run_diamond_screen(query_fa, args.dmnd_db, hits_tsv, args.threads)

    summary = {
        "viral_families_genera_flagged": len(rank_rows),
        "cdna_unmapped_k2u_total": cdna_tso + len(cdna_kept),
        "cdna_tso_primer_filtered": cdna_tso,
        "cdna_unknown_screened": len(cdna_kept),
        "dna_unknown_screened": len(dna_kept),
        "total_unknown_screened": len(cdna_kept) + len(dna_kept),
        "diamond_any_hit_reads": raw_q_hits,
        "diamond_strong_hit_reads_aa50_bits50": len(strong_hits),
        "strong_hit_titles": dict(collections.Counter(h["stitle"] for h in strong_hits).most_common(10)),
    }
    sum_json = os.path.join(args.outdir, "unknown_screen_summary.json")
    with open(sum_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

