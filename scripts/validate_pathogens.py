#!/usr/bin/env python3
"""Competitive reference-based pathogen validation for non-host Nanopore reads.

Builds a combined competitive FASTA (human masked reference + candidate pathogen
references with prefixed contig names), maps DNA and cDNA non-host libraries
with minimap2, filters out unmapped/human primary alignments into sorted BAMs,
and summarizes per-species and per-contig coverage, depth, MAPQ, and identity.
"""

import argparse
import collections
import gzip
import json
import os
import subprocess
import sys


SPECIES_REFS = [
    ("FLU_H3N2", "influenza_a_h3n2_ny392.fna.gz"),
    ("FLU_H1N1", "influenza_a_h1n1_ca07.fna.gz"),
    ("CANDIDA_ALBICANS", "candida_albicans_sc5314.fna.gz"),
    ("MALASSEZIA_RESTRICTA", "malassezia_restricta.fna.gz"),
    ("ECOLI", "escherichia_coli_mg1655.fna.gz"),
    ("KLEBSIELLA_PNEUMONIAE", "klebsiella_pneumoniae_hs11286.fna.gz"),
    ("ACINETOBACTER_PITTII", "acinetobacter_pittii.fna.gz"),
    ("VIBRIO_ALGINOLYTICUS", "vibrio_alginolyticus.fna.gz"),
]


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def build_combined_reference(human_ref, ref_dir, out_fa_gz):
    """Concatenate human and pathogen references with `<TAG>|` contig prefixes."""
    if os.path.exists(out_fa_gz) and os.path.getsize(out_fa_gz) > 0:
        print(f"[ref] Reusing existing {out_fa_gz}", flush=True)
        return

    tmp_fa = out_fa_gz[:-3] if out_fa_gz.endswith(".gz") else out_fa_gz
    sources = [("HUMAN", human_ref)] + [
        (tag, os.path.join(ref_dir, fname)) for tag, fname in SPECIES_REFS
    ]
    with open(tmp_fa, "wt") as out:
        for tag, path in sources:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Missing reference: {path}")
            print(f"[ref] Adding {tag} from {path}", flush=True)
            with open_maybe_gz(path, "rt") as inp:
                for line in inp:
                    if line.startswith(">"):
                        header = line[1:].strip()
                        contig_id = header.split()[0]
                        rest = header[len(contig_id):]
                        out.write(f">{tag}|{contig_id}{rest}\n")
                    else:
                        out.write(line)

    print(f"[ref] Compressing {tmp_fa} -> {out_fa_gz}", flush=True)
    subprocess.run(["bgzip", "-f", "-@", "8", tmp_fa], check=True)


def build_minimap2_index(ref_fa_gz, out_mmi, threads):
    if os.path.exists(out_mmi) and os.path.getsize(out_mmi) > 0:
        print(f"[index] Reusing existing {out_mmi}", flush=True)
        return
    print(f"[index] Building {out_mmi}", flush=True)
    subprocess.run(
        ["minimap2", "-x", "map-ont", "-t", str(threads), "-d", out_mmi, ref_fa_gz],
        check=True,
    )


def parse_cigar_aligned_len(cigar):
    """Return reference-aligned length (M, D, N, =, X) and total query aligned (M, I, =, X)."""
    ref_len = 0
    qry_len = 0
    num = 0
    for ch in cigar:
        if "0" <= ch <= "9":
            num = num * 10 + (ord(ch) - 48)
        else:
            if ch in ("M", "=", "X"):
                ref_len += num
                qry_len += num
            elif ch in ("D", "N"):
                ref_len += num
            elif ch == "I":
                qry_len += num
            num = 0
    return ref_len, qry_len


def map_and_filter(fastq_gz, mmi, out_bam, threads):
    """Map reads to combined competitive index and keep non-human primary alignments."""
    if os.path.exists(out_bam) and os.path.exists(out_bam + ".bai"):
        stats_json = out_bam + ".map_stats.json"
        if os.path.exists(stats_json):
            print(f"[map] Reusing existing {out_bam}", flush=True)
            with open(stats_json) as f:
                return json.load(f)

    print(f"[map] Mapping {fastq_gz} -> {out_bam}", flush=True)
    mm2_cmd = [
        "minimap2",
        "-ax",
        "map-ont",
        "--secondary=no",
        "-t",
        str(threads),
        mmi,
        fastq_gz,
    ]
    unsorted_bam = out_bam + ".unsorted.bam"
    mm2_proc = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, text=True)
    sam_view_proc = subprocess.Popen(
        ["samtools", "view", "-b", "-@", "4", "-o", unsorted_bam, "-"],
        stdin=subprocess.PIPE,
        text=True,
    )

    counts = {
        "total_primary": 0,
        "unmapped": 0,
        "human_primary": 0,
        "nonhuman_primary": 0,
    }

    for line in mm2_proc.stdout:
        if line.startswith("@"):
            if line.startswith("@SQ") and "\tSN:HUMAN|" in line:
                continue
            sam_view_proc.stdin.write(line)
            continue

        fields = line.split("\t", 6)
        flag = int(fields[1])
        if flag & 0x900:
            # skip secondary/supplementary if any
            continue
        counts["total_primary"] += 1
        if flag & 0x4:
            counts["unmapped"] += 1
            continue
        rname = fields[2]
        if rname.startswith("HUMAN|"):
            counts["human_primary"] += 1
            continue

        counts["nonhuman_primary"] += 1
        sam_view_proc.stdin.write(line)

    sam_view_proc.stdin.close()
    mm2_ret = mm2_proc.wait()
    sv_ret = sam_view_proc.wait()
    if mm2_ret != 0 or sv_ret != 0:
        raise RuntimeError(f"minimap2 ({mm2_ret}) or samtools view ({sv_ret}) failed")

    subprocess.run(
        ["samtools", "sort", "-@", "8", "-o", out_bam, unsorted_bam], check=True
    )
    subprocess.run(["samtools", "index", out_bam], check=True)
    os.remove(unsorted_bam)

    with open(out_bam + ".map_stats.json", "w") as f:
        json.dump(counts, f, indent=2)
    return counts


def summarize_bam(bam_path, lib_name):
    """Compute per-species and per-contig read + coverage stats from filtered BAM."""
    # Reference lengths from BAM header
    hdr = subprocess.check_output(["samtools", "view", "-H", bam_path], text=True)
    contig_len = {}
    for line in hdr.splitlines():
        if line.startswith("@SQ"):
            sn = ln = None
            for part in line.split("\t")[1:]:
                if part.startswith("SN:"):
                    sn = part[3:]
                elif part.startswith("LN:"):
                    ln = int(part[3:])
            if sn and ln is not None:
                contig_len[sn] = ln

    # Read stats via samtools view
    sp_stats = collections.defaultdict(
        lambda: {
            "reads": 0,
            "reads_mapq20": 0,
            "aligned_bases": 0,
            "nm_sum": 0,
            "mapq_sum": 0,
        }
    )
    contig_reads = collections.Counter()
    contig_reads_q20 = collections.Counter()

    view_proc = subprocess.Popen(
        ["samtools", "view", bam_path], stdout=subprocess.PIPE, text=True
    )
    for line in view_proc.stdout:
        parts = line.rstrip("\n").split("\t")
        rname = parts[2]
        mapq = int(parts[4])
        cigar = parts[5]
        species = rname.split("|", 1)[0]

        _, qry_aln = parse_cigar_aligned_len(cigar)
        nm = 0
        for tag in parts[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag[5:])
                break

        st = sp_stats[species]
        st["reads"] += 1
        if mapq >= 20:
            st["reads_mapq20"] += 1
            contig_reads_q20[rname] += 1
        st["aligned_bases"] += qry_aln
        st["nm_sum"] += nm
        st["mapq_sum"] += mapq
        contig_reads[rname] += 1

    if view_proc.wait() != 0:
        raise RuntimeError(f"samtools view failed on {bam_path}")

    # Coverage stats via samtools coverage (fast, per-contig)
    cov_out = subprocess.check_output(["samtools", "coverage", bam_path], text=True)
    sp_cov = collections.defaultdict(
        lambda: {"ref_len": 0, "covered_bases": 0, "depth_weighted_sum": 0.0, "max_cov_pct": 0.0}
    )
    contig_rows = []
    for line in cov_out.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        rname = cols[0]
        endpos = int(cols[2])
        numreads = int(cols[3])
        covbases = int(cols[4])
        coverage_pct = float(cols[5])
        meandepth = float(cols[6])
        meanmapq = float(cols[8])
        species = rname.split("|", 1)[0]

        sp_cov[species]["ref_len"] += endpos
        sp_cov[species]["covered_bases"] += covbases
        sp_cov[species]["depth_weighted_sum"] += meandepth * endpos
        if coverage_pct > sp_cov[species]["max_cov_pct"]:
            sp_cov[species]["max_cov_pct"] = coverage_pct

        if numreads > 0:
            contig_rows.append(
                {
                    "library": lib_name,
                    "species": species,
                    "contig": rname.split("|", 1)[1],
                    "length": endpos,
                    "reads": numreads,
                    "reads_mapq20": contig_reads_q20[rname],
                    "covered_bases": covbases,
                    "coverage_pct": round(coverage_pct, 4),
                    "mean_depth": round(meandepth, 4),
                    "mean_mapq": round(meanmapq, 2),
                }
            )

    species_rows = []
    for species, _ in SPECIES_REFS:
        st = sp_stats[species]
        cv = sp_cov[species]
        ref_len = cv["ref_len"] or sum(
            ln for sn, ln in contig_len.items() if sn.startswith(species + "|")
        )
        reads = st["reads"]
        aln_bases = st["aligned_bases"]
        mean_ident = (1.0 - st["nm_sum"] / aln_bases) * 100.0 if aln_bases > 0 else 0.0
        mean_mapq = st["mapq_sum"] / reads if reads > 0 else 0.0
        breadth = (cv["covered_bases"] / ref_len * 100.0) if ref_len > 0 else 0.0
        mean_depth = (cv["depth_weighted_sum"] / ref_len) if ref_len > 0 else 0.0
        species_rows.append(
            {
                "library": lib_name,
                "species": species,
                "ref_len": ref_len,
                "mapped_reads": reads,
                "mapped_reads_mapq20": st["reads_mapq20"],
                "mean_mapq": round(mean_mapq, 2),
                "mean_identity_pct": round(mean_ident, 2),
                "covered_bases": cv["covered_bases"],
                "breadth_pct": round(breadth, 4),
                "mean_depth": round(mean_depth, 4),
                "max_contig_breadth_pct": round(cv["max_cov_pct"], 2),
            }
        )

    return species_rows, contig_rows


def write_tsv(path, rows):
    if not rows:
        return
    cols = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--outdir", default="analysis/pathogen_validation")
    args = p.parse_args()

    ref_dir = os.path.join(args.outdir, "references")
    combined_fa_gz = os.path.join(ref_dir, "combined_competitive.fa.gz")
    mmi = os.path.join(ref_dir, "combined_competitive.mmi")

    build_combined_reference("references/human-masked.fa.gz", ref_dir, combined_fa_gz)
    build_minimap2_index(combined_fa_gz, mmi, args.threads)

    libs = [
        ("DNA_SG26-GS01", "clean_results/sample_001.SG26-GS01.nonhost.fastq.gz"),
        ("cDNA_ST26-UM01", "clean_results/sample_001.ST26-UM01.nonhost.fastq.gz"),
    ]

    all_map_stats = {}
    all_species_rows = []
    all_contig_rows = []

    for lib_name, fq in libs:
        bam = os.path.join(args.outdir, f"{lib_name}.nonhuman.bam")
        mstats = map_and_filter(fq, mmi, bam, args.threads)
        all_map_stats[lib_name] = mstats
        sp_rows, ctg_rows = summarize_bam(bam, lib_name)
        all_species_rows.extend(sp_rows)
        all_contig_rows.extend(ctg_rows)

    with open(os.path.join(args.outdir, "mapping_overview.json"), "w") as f:
        json.dump(all_map_stats, f, indent=2)
    write_tsv(os.path.join(args.outdir, "species_validation_summary.tsv"), all_species_rows)
    write_tsv(os.path.join(args.outdir, "contig_validation_summary.tsv"), all_contig_rows)
    print("[done] Validation summaries written to", args.outdir, flush=True)


if __name__ == "__main__":
    main()

