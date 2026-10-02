#!/usr/bin/env python3
"""
10_gcf_prevalence.py

Purpose: Compute GCF prevalence across the 110 genomes for each BiG-SCAPE cutoff.

Inputs:  results/07_bigscape/raw_output/  (BiG-SCAPE network files)
         results/04_gtdbtk/final_genome_list.txt
Outputs:
  - results/07_bigscape/parsed/gcf_prevalence_c0.3.tsv
  - results/07_bigscape/parsed/gcf_prevalence_c0.5.tsv
  - results/07_bigscape/parsed/gcf_prevalence_c0.7.tsv
  - results/07_bigscape/parsed/gcf_summary.txt
  - results/07_bigscape/parsed/gcf_overlap_devosia_A.tsv
"""
from pathlib import Path
from collections import defaultdict
import re

ROOT = Path(__file__).resolve().parent.parent
BS = ROOT / "results" / "07_bigscape"
PARSED = BS / "parsed"
GTDBTK = ROOT / "results" / "04_gtdbtk"

# Prevalence thresholds (matching manuscript)
CORE_THRESHOLD = 90.0
RARE_THRESHOLD = 10.0


def parse_gcf_clustering(network_file, cutoff):
    """
    Parse BiG-SCAPE network file to get GCF assignments.
    BiG-SCAPE 2.0 output format: TSV with BGC and family columns.
    """
    assignments = {}  # bgc_id -> gcf_id
    with open(network_file) as f:
        header = f.readline().strip().split("\t")
        # Typical columns: BGC, GCF, Family, ...
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            bgc_id = parts[0]
            # BiG-SCAPE assigns cluster family names like "BGC0001_1"
            gcf_id = parts[1] if len(parts) > 1 else None
            if gcf_id and gcf_id != "NULL":
                assignments[bgc_id] = gcf_id
    return assignments


def bgc_to_accession(bgc_id):
    """Extract genome accession from bgc_id like GCF_000001_region001."""
    # bgc_id format: <accession>_region<NNN>
    m = re.match(r"^(.+?)_region\d+$", bgc_id)
    if m:
        return m.group(1)
    return bgc_id


def main():
    print("=" * 60)
    print("10_gcf_prevalence.py")
    print("=" * 60)

    # Load final genome list
    final_list = GTDBTK / "final_genome_list.txt"
    if not final_list.exists():
        print(f"ERROR: {final_list} not found")
        return
    genomes = [g.strip() for g in final_list.read_text().splitlines() if g.strip()]
    n_genomes = len(genomes)
    print(f"Loaded {n_genomes} final genomes")

    # Load Devosia_A genomes (from GTDB-Tk output)
    devosia_a = set()
    tax_file = GTDBTK / "taxonomy_assignments.tsv"
    if tax_file.exists():
        with open(tax_file) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[2] == "Devosia_A":
                    devosia_a.add(parts[0])
    print(f"Devosia_A genomes: {len(devosia_a)}")

    for cutoff in ["0.3", "0.5", "0.7"]:
        print(f"\n--- Cutoff c{cutoff} ---")

        # Find BiG-SCAPE output network files for this cutoff
        # BiG-SCAPE output: raw_output/<bin>/<cutoff>/network_files/*_clustering_<cutoff>.tsv
        network_files = list(BS.rglob(f"*_clustering_c{cutoff}.tsv"))
        if not network_files:
            # Alternative: use the mix directory
            network_files = list(BS.rglob(f"*c{cutoff}*.tsv"))
            network_files = [f for f in network_files if "network" in str(f).lower()]

        print(f"  Found {len(network_files)} network files")

        assignments = {}
        for nf in network_files:
            assignments.update(parse_gcf_clustering(nf, cutoff))

        print(f"  Total BGC assignments: {len(assignments)}")

        # Build GCF -> genome set
        gcf_genomes = defaultdict(set)
        gcf_classes = defaultdict(set)

        for bgc_id, gcf_id in assignments.items():
            acc = bgc_to_accession(bgc_id)
            gcf_genomes[gcf_id].add(acc)

        # Compute prevalence and classify
        rows = []
        for gcf_id, accs in gcf_genomes.items():
            # Only count genomes in our final list
            accs_in_list = {a for a in accs if a in genomes}
            n_present = len(accs_in_list)
            prevalence = 100 * n_present / n_genomes

            if prevalence >= CORE_THRESHOLD:
                category = "core"
            elif prevalence < RARE_THRESHOLD:
                category = "rare"
            else:
                category = "intermediate"

            rows.append({
                "gcf_id": gcf_id,
                "n_genomes": n_present,
                "prevalence": round(prevalence, 2),
                "category": category,
                "genomes": ",".join(sorted(accs_in_list)),
            })

        # Sort by prevalence descending
        rows.sort(key=lambda r: -r["prevalence"])

        # Write table
        out_file = PARSED / f"gcf_prevalence_c{cutoff}.tsv"
        with open(out_file, "w") as f:
            f.write("gcf_id\tn_genomes\tprevalence\tcategory\tgenomes\n")
            for r in rows:
                f.write(f"{r['gcf_id']}\t{r['n_genomes']}\t{r['prevalence']}\t{r['category']}\t{r['genomes']}\n")

        # Summary
        n_core = sum(1 for r in rows if r["category"] == "core")
        n_inter = sum(1 for r in rows if r["category"] == "intermediate")
        n_rare = sum(1 for r in rows if r["category"] == "rare")
        print(f"  GCFs: {len(rows)}  |  core={n_core}, intermediate={n_inter}, rare={n_rare}")

    # Devosia_A overlap analysis at c0.7
    print("\n--- Devosia_A overlap at c0.7 ---")
    c07_file = PARSED / "gcf_prevalence_c0.7.tsv"
    if c07_file.exists() and devosia_a:
        dev_gcfs = set()
        da_gcfs = set()
        shared = set()
        with open(c07_file) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split("\t")
                gcf_id = parts[0]
                accs = set(parts[4].split(","))
                in_dev = bool(accs - devosia_a)
                in_da = bool(accs & devosia_a)
                if in_dev and in_da:
                    shared.add(gcf_id)
                elif in_dev:
                    dev_gcfs.add(gcf_id)
                elif in_da:
                    da_gcfs.add(gcf_id)

        overlap_file = PARSED / "gcf_overlap_devosia_A.tsv"
        with open(overlap_file, "w") as f:
            f.write("category\tn_gcfs\n")
            f.write(f"Devosia_specific\t{len(dev_gcfs)}\n")
            f.write(f"Devosia_A_specific\t{len(da_gcfs)}\n")
            f.write(f"shared\t{len(shared)}\n")
            f.write(f"total\t{len(dev_gcfs) + len(da_gcfs) + len(shared)}\n")

        print(f"  Devosia-specific:   {len(dev_gcfs)}")
        print(f"  Devosia_A-specific: {len(da_gcfs)}")
        print(f"  Shared:             {len(shared)}")
        print(f"  Total:              {len(dev_gcfs) + len(da_gcfs) + len(shared)}")


if __name__ == "__main__":
    main()
