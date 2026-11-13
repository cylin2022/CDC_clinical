#!/usr/bin/env python3
"""Collect Kraken2 reports (with --report-minimizer-data) into species/genus tables."""
import csv
import sys
from pathlib import Path

RANKS = {'S': 'species', 'G': 'genus'}
MIN_READS_REPLICATE = 3
# Genera reported as reagent/kit contaminants (Salter et al. 2014 BMC Biol 12:87, Table 1;
# Eisenhofer et al. 2019 Trends Microbiol 27:105). A flag means "interpret with caution",
# NOT exclusion: several (e.g. Streptococcus, Escherichia, Pseudomonas) are also true pathogens.
CONTAMINANT_GENERA = set("""
Afipia Aquabacterium Asticcacaulis Aurantimonas Beijerinckia Bosea Bradyrhizobium Brevundimonas
Caulobacter Craurococcus Devosia Hoeflea Mesorhizobium Methylobacterium Novosphingobium
Ochrobactrum Paracoccus Pedomicrobium Phyllobacterium Rhizobium Roseomonas Sphingobium
Sphingomonas Sphingopyxis Acidovorax Azoarcus Burkholderia Comamonas Cupriavidus Curvibacter
Delftia Duganella Herbaspirillum Janthinobacterium Kingella Leptothrix Limnobacter Massilia
Methylophilus Methyloversatilis Oxalobacter Pelomonas Polaromonas Ralstonia Schlegelella
Sulfuritalea Undibacterium Variovorax Acinetobacter Enhydrobacter Enterobacter Escherichia
Nevskia Pseudomonas Pseudoxanthomonas Psychrobacter Stenotrophomonas Xanthomonas Aeromicrobium
Arthrobacter Beutenbergia Brevibacterium Corynebacterium Curtobacterium Dietzia
Geodermatophilus Janibacter Kocuria Microbacterium Micrococcus Microlunatus Patulibacter
Propionibacterium Cutibacterium Rhodococcus Tsukamurella Abiotrophia Bacillus Brevibacillus
Brochothrix Facklamia Paenibacillus Streptococcus Chryseobacterium Dyadobacter Flavobacterium
Hydrotalea Niastella Olivibacter Pedobacter Wautersiella Deinococcus
""".split())


def source_libraries(out):
    """Map report name 'source.<source_id>' -> library_type using source_summary.tsv."""
    path = out.parent / 'source_summary.tsv'
    if not path.is_file():
        return {}
    with path.open() as fh:
        return {f"source.{r['source_id']}": r['library_type']
                for r in csv.DictReader(fh, delimiter='\t')}


def parse(path):
    rows, total = [], 0
    for line in path.read_text().splitlines():
        f = line.split('\t')
        if len(f) < 8:
            continue
        pct, clade, direct, minimizers, distinct, rank, taxid, name = f[:8]
        clade = int(clade)
        if rank in ('U', 'R'):
            total += clade
        rows.append(dict(rank=rank, taxid=taxid, name=name.strip(), clade_reads=clade,
                         direct_reads=int(direct), minimizers=int(minimizers),
                         distinct_minimizers=int(distinct)))
    return rows, total


def main(out_dir):
    out = Path(out_dir)
    reports = sorted(out.glob('*.report.txt'))
    if not reports:
        raise SystemExit(f'No reports in {out}')
    overview = []
    libs = source_libraries(out)
    lib_types = sorted(set(libs.values()))
    for code, label in RANKS.items():
        table = {}
        for rep in reports:
            name = rep.name[:-len('.report.txt')]
            rows, total = parse(rep)
            if code == 'S':
                unclassified = next((r['clade_reads'] for r in rows if r['rank'] == 'U'), 0)
                overview.append(dict(report=name, total_reads=total,
                                     unclassified=unclassified,
                                     classified=total - unclassified))
            for r in rows:
                if r['rank'] != code:
                    continue
                key = (r['taxid'], r['name'])
                table.setdefault(key, {})[name] = (r['clade_reads'], r['distinct_minimizers'], total)
        names = [rep.name[:-len('.report.txt')] for rep in reports]
        path = out / f'{label}_table.tsv'
        with path.open('w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            header = ['taxid', 'name']
            header += ['contaminant_genus_flag'] + [f'both_runs_{t}_ge{MIN_READS_REPLICATE}reads'
                                                    for t in lib_types]
            for n in names:
                header += [f'{n}.reads', f'{n}.pct_of_nonhost', f'{n}.distinct_minimizers']
            w.writerow(header)
            order = sorted(table.items(), key=lambda kv: -max(v[0] for v in kv[1].values()))
            for (taxid, sci), per in order:
                genus = sci.split()[0] if sci else ''
                line = [taxid, sci, 'yes' if genus in CONTAMINANT_GENERA else '']
                for t in lib_types:
                    runs = [per.get(n, (0,))[0] for n, lt in libs.items() if lt == t]
                    line.append('yes' if runs and all(x >= MIN_READS_REPLICATE for x in runs) else 'no')
                for n in names:
                    reads, distinct, total = per.get(n, (0, 0, 0))
                    line += [reads, f'{100 * reads / total:.4f}' if total else '0', distinct]
                w.writerow(line)
    with (out / 'overview.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(overview[0]), delimiter='\t')
        w.writeheader()
        w.writerows(overview)
    print(f'Wrote overview.tsv, species_table.tsv, genus_table.tsv in {out}')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'clean_results/kraken2')

