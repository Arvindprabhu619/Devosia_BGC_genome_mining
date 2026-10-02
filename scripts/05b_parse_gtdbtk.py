#!/usr/bin/env python3
"""05b_parse_gtdbtk.py — Parse GTDB-Tk summary into taxonomy tables."""
import csv
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
GTDBTK_OUT = ROOT / "results" / "04_gtdbtk" / "gtdbtk_output" / "classify"
RES = ROOT / "results" / "04_gtdbtk"

summary = GTDBTK_OUT / "gtdbtk.bac120.summary.tsv"

rows = []
with open(summary) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)

print(f"Total taxonomy rows: {len(rows)}")
print()

# Extract genus from taxonomy string
def get_genus(tax):
    for part in tax.split(";"):
        if part.startswith("g__"):
            return part.replace("g__", "")
    return "Unknown"

genus_counts = Counter()
devosia, devosia_a, other = [], [], []

for row in rows:
    acc = row["user_genome"].replace(".fna", "").replace("query_", "")
    tax = row["classification"]
    genus = get_genus(tax)
    genus_counts[genus] += 1

    if genus == "Devosia_A":
        devosia_a.append(acc)
    elif genus == "Devosia":
        devosia.append(acc)
    else:
        other.append((acc, genus))

print("=== Genus distribution across all 133 genomes ===")
for g, c in genus_counts.most_common():
    print(f"  {g:30s} {c}")

# Write outputs
with open(RES / "taxonomy_assignments.tsv", "w") as f:
    f.write("accession\ttaxonomy\tgenus\tstatus\n")
    for row in rows:
        acc = row["user_genome"].replace(".fna", "").replace("query_", "")
        tax = row["classification"]
        genus = get_genus(tax)
        if genus == "Devosia":
            status = "Devosia"
        elif genus == "Devosia_A":
            status = "Devosia_A"
        else:
            status = "Excluded"
        f.write(f"{acc}\t{tax}\t{genus}\t{status}\n")

with open(RES / "final_genome_list.txt", "w") as f:
    for acc in devosia + devosia_a:
        f.write(acc + "\n")

with open(RES / "excluded_taxa.txt", "w") as f:
    for acc, genus in other:
        f.write(f"{acc}\t{genus}\n")

print()
print("=== Final classification ===")
print(f"Devosia:       {len(devosia)}")
print(f"Devosia_A:     {len(devosia_a)}")
print(f"Excluded:      {len(other)}")
print(f"FINAL:         {len(devosia) + len(devosia_a)}")
