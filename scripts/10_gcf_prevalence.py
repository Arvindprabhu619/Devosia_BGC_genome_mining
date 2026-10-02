#!/usr/bin/env python3
# ============================================================
# scripts/10_gcf_prevalence.py
# ============================================================
# Purpose: Compute GCF prevalence across genomes for each cutoff.
#
# Inputs:   results/09_bigscape/raw_output/
#           results/05_gtdbtk/final_genome_list.txt
# Outputs:
#   results/09_bigscape/parsed/gcf_prevalence_c{0.3,0.5,0.7}.tsv
#   results/09_bigscape/parsed/gcf_overlap_devosia_A.tsv
#
# Thresholds:
#   Core:         prevalence >= 90%
#   Intermediate: 10% <= prevalence < 90%
#   Rare:         prevalence < 10%
#
# Runtime:  15 min
# ============================================================
from pathlib import Path
from collections import defaultdict
import re
import csv

ROOT = Path(__file__).resolve().parent.parent
BS = ROOT / "results" / "09_bigscape"
PARSED = BS / "parsed"
GTDBTK = ROOT / "results" / "05_gtdbtk"

CORE_THRESHOLD = 90.0
RARE_THRESHOLD = 10.0


def bgc_to_acc(bgc_id):
    m = re.match(r"^(.+?)_region\d+$", bgc_id)
    return m.group(1) if m else bgc_id


def parse_cluster_file(f):
    """Parse a BiG-SCAPE cluster file (TSV) into {bgc_id: gcf_id}."""
    out = {}
    with open(f) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            out[parts[0]] = parts[1]
    return out


def main():
    print("=" * 60)
    print("10_gcf_prevalence.py")
    print("=" * 60)

    final_list = GTDBTK / "final_genome_list.txt"
    genomes = [g.strip() for g in final_list.read_text().splitlines() if g.strip()]
    n_genomes = len(genomes)
    print(f"Genomes: {n_genomes}")

    devosia_a = set()
    tax_file = GTDBTK / "taxonomy_assignments.tsv"
    if tax_file.exists():
        with open(tax_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("status") == "Devosia_A":
                    devosia_a.add(row["accession"])
    print(f"Devosia_A: {len(devosia_a)}")

    for cutoff in ["0.3", "0.5", "0.7"]:
        print(f"\n--- Cutoff c{cutoff} ---")

        # Find cluster TSVs from BiG-SCAPE output
        cluster_files = list(BS.rglob(f"*_clustering_c{cutoff}.tsv"))

        if not cluster_files:
            print(f"  No cluster file found for cutoff {cutoff}")
            continue

        assignments = {}
        for cf in cluster_files:
            assignments.update(parse_cluster_file(cf))

        print(f"  BGC assignments: {len(assignments)}")

        # Map GCF -> genomes
        gcf_genomes = defaultdict(set)
        for bgc_id, gcf_id in assignments.items():
            acc = bgc_to_acc(bgc_id)
            gcf_genomes[gcf_id].add(acc)

        rows = []
        for gcf_id, accs in gcf_genomes.items():
            accs_final = {a for a in accs if a in genomes}
            n = len(accs_final)
            prev = 100.0 * n / n_genomes
            if prev >= CORE_THRESHOLD:
                cat = "core"
            elif prev < RARE_THRESHOLD:
                cat = "rare"
            else:
                cat = "intermediate"
            rows.append({
                "gcf_id": gcf_id,
                "n_genomes": n,
                "prevalence": round(prev, 2),
                "category": cat,
                "genomes": ",".join(sorted(accs_final)),
            })

        rows.sort(key=lambda r: -r["prevalence"])

        with open(PARSED / f"gcf_prevalence_c{cutoff}.tsv", "w") as f:
            f.write("gcf_id\tn_genomes\tprevalence\tcategory\tgenomes\n")
            for r in rows:
                f.write(f"{r['gcf_id']}\t{r['n_genomes']}\t{r['prevalence']}\t{r['category']}\t{r['genomes']}\n")

        nc = sum(1 for r in rows if r["category"] == "core")
        ni = sum(1 for r in rows if r["category"] == "intermediate")
        nr = sum(1 for r in rows if r["category"] == "rare")
        print(f"  GCFs: {len(rows)} | core={nc}, intermediate={ni}, rare={nr}")

    # Devosia vs Devosia_A overlap at c0.7
    c07 = PARSED / "gcf_prevalence_c0.7.tsv"
    if c07.exists() and devosia_a:
        dev_only, da_only, shared = set(), set(), set()
        with open(c07) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                accs = set(row["genomes"].split(",")) if row["genomes"] else set()
                in_dev = bool(accs - devosia_a)
                in_da = bool(accs & devosia_a)
                if in_dev and in_da:
                    shared.add(row["gcf_id"])
                elif in_dev:
                    dev_only.add(row["gcf_id"])
                elif in_da:
                    da_only.add(row["gcf_id"])

        with open(PARSED / "gcf_overlap_devosia_A.tsv", "w") as f:
            f.write("category\tn_gcfs\n")
            f.write(f"Devosia_specific\t{len(dev_only)}\n")
            f.write(f"Devosia_A_specific\t{len(da_only)}\n")
            f.write(f"shared\t{len(shared)}\n")
            f.write(f"total\t{len(dev_only) + len(da_only) + len(shared)}\n")

        print(f"\n--- Devosia vs Devosia_A at c0.7 ---")
        print(f"  Devosia-specific:   {len(dev_only)}")
        print(f"  Devosia_A-specific: {len(da_only)}")
        print(f"  Shared:             {len(shared)}")


if __name__ == "__main__":
    main()
