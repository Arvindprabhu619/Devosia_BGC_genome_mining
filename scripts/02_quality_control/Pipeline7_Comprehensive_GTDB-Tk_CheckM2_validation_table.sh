
conda activate gtdbtk

export GTDBTK_DATA_PATH="$HOME/Microbial_minning_project/gtdbtk_db/release232"

BASE="$PROJECT_ROOT"
QC="$BASE/Devosia_CheckM2_QC.tsv"
GTDB="$BASE/gtdbtk_r232/classify/gtdbtk.bac120.summary.tsv"
OUT="$BASE/validation"

mkdir -p "$OUT"

echo "============================================================"
echo "DEVOSIA GENOME VALIDATION"
echo "CheckM2 + GTDB-Tk R232"
echo "============================================================"

echo
echo "Date:"
date

echo
echo "=== Input files ==="
ls -lh "$QC"
ls -lh "$GTDB"

echo
echo "=== Genome count ==="
find "$BASE/genomes" -maxdepth 1 -type f -name "*.fna" | wc -l

echo
echo "=== CheckM2 rows ==="
wc -l "$QC"

echo
echo "=== GTDB-Tk rows ==="
wc -l "$GTDB"

echo
echo "============================================================"
echo "GENERATING INTEGRATED VALIDATION TABLE"
echo "============================================================"

python - "$QC" "$GTDB" "$OUT" <<'PY'

import csv
import sys
from pathlib import Path
from collections import Counter

qc_file = Path(sys.argv[1])
gtdb_file = Path(sys.argv[2])
out_dir = Path(sys.argv[3])

out_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Read CheckM2
# ------------------------------------------------------------

with open(qc_file, newline="") as fh:
    qc_rows = list(csv.DictReader(fh, delimiter="\t"))

qc = {}

for r in qc_rows:
    name = r["Name"].strip()
    qc[name] = r

# ------------------------------------------------------------
# Read GTDB-Tk
# ------------------------------------------------------------

with open(gtdb_file, newline="") as fh:
    gtdb_rows = list(csv.DictReader(fh, delimiter="\t"))

# ------------------------------------------------------------
# Taxonomy parser
# ------------------------------------------------------------

def parse_taxonomy(tax):

    result = {
        "domain": "",
        "phylum": "",
        "class": "",
        "order": "",
        "family": "",
        "genus": "",
        "species": ""
    }

    if not tax:
        return result

    for item in tax.split(";"):

        item = item.strip()

        if item.startswith("d__"):
            result["domain"] = item[3:]

        elif item.startswith("p__"):
            result["phylum"] = item[3:]

        elif item.startswith("c__"):
            result["class"] = item[3:]

        elif item.startswith("o__"):
            result["order"] = item[3:]

        elif item.startswith("f__"):
            result["family"] = item[3:]

        elif item.startswith("g__"):
            result["genus"] = item[3:]

        elif item.startswith("s__"):
            result["species"] = item[3:]

    return result


# ------------------------------------------------------------
# Build integrated table
# ------------------------------------------------------------

integrated = []

for r in gtdb_rows:

    genome = r["user_genome"].strip()

    # GTDB genome names do not necessarily contain the exact
    # CheckM2 suffix. Match using accession where possible.
    qc_match = None

    if genome in qc:
        qc_match = qc[genome]

    else:
        accession = genome.split("_")[0] + "_" + genome.split("_")[1] \
            if genome.startswith(("GCA_", "GCF_")) else genome

        candidates = [
            q for name, q in qc.items()
            if name.startswith(accession)
        ]

        if len(candidates) == 1:
            qc_match = candidates[0]

    tax = parse_taxonomy(r["classification"])

    row = {
        "Genome": genome,

        "CheckM2_Completeness":
            qc_match["Completeness"] if qc_match else "",

        "CheckM2_Contamination":
            qc_match["Contamination"] if qc_match else "",

        "CheckM2_QC_status":
            qc_match["QC_status"] if qc_match else "",

        "GTDB_Domain": tax["domain"],
        "GTDB_Phylum": tax["phylum"],
        "GTDB_Class": tax["class"],
        "GTDB_Order": tax["order"],
        "GTDB_Family": tax["family"],
        "GTDB_Genus": tax["genus"],
        "GTDB_Species": tax["species"],

        "GTDB_Classification_Method":
            r["classification_method"],

        "GTDB_Note":
            r["note"],

        "Closest_Genome":
            r["closest_genome_reference"],

        "Closest_Genome_ANI":
            r["closest_genome_ani"],

        "Closest_Genome_AF":
            r["closest_genome_af"],

        "MSA_Percent":
            r["msa_percent"],

        "Translation_Table":
            r["translation_table"],

        "RED_Value":
            r["red_value"],

        "Warnings":
            r["warnings"],
    }

    integrated.append(row)


# ------------------------------------------------------------
# Write integrated table
# ------------------------------------------------------------

fields = list(integrated[0].keys())

integrated_file = out_dir / "Devosia_GTDB_CheckM2_integrated.tsv"

with open(integrated_file, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(integrated)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("=" * 70)
print("INTEGRATED VALIDATION SUMMARY")
print("=" * 70)

print(f"CheckM2 genomes       : {len(qc_rows)}")
print(f"GTDB-Tk genomes       : {len(gtdb_rows)}")
print(f"Integrated genomes    : {len(integrated)}")

# ------------------------------------------------------------
# CheckM2 status
# ------------------------------------------------------------

print()
print("=== CheckM2 QC STATUS ===")

qc_status = Counter(
    r["CheckM2_QC_status"] or "MISSING"
    for r in integrated
)

for k, v in qc_status.most_common():
    print(f"{k:20s} {v}")

# ------------------------------------------------------------
# GTDB genus
# ------------------------------------------------------------

print()
print("=== GTDB-Tk GENUS DISTRIBUTION ===")

genera = Counter(
    r["GTDB_Genus"] or "UNCLASSIFIED"
    for r in integrated
)

for k, v in genera.most_common():
    print(f"{k:30s} {v}")

# ------------------------------------------------------------
# GTDB family
# ------------------------------------------------------------

print()
print("=== GTDB-Tk FAMILY DISTRIBUTION ===")

families = Counter(
    r["GTDB_Family"] or "UNCLASSIFIED"
    for r in integrated
)

for k, v in families.most_common():
    print(f"{k:30s} {v}")

# ------------------------------------------------------------
# Classification methods
# ------------------------------------------------------------

print()
print("=== CLASSIFICATION METHODS ===")

methods = Counter(
    r["GTDB_Classification_Method"] or "MISSING"
    for r in integrated
)

for k, v in methods.most_common():
    print(f"{k:50s} {v}")

# ------------------------------------------------------------
# Species assignments
# ------------------------------------------------------------

species_assigned = sum(
    bool(r["GTDB_Species"]) and r["GTDB_Species"] != ""
    for r in integrated
)

print()
print("=== SPECIES ASSIGNMENT ===")
print(f"Species assigned     : {species_assigned}")
print(f"Species unassigned   : {len(integrated) - species_assigned}")

# ------------------------------------------------------------
# Devosia
# ------------------------------------------------------------

devosia = [
    r for r in integrated
    if r["GTDB_Genus"].lower().startswith("devosia")
]

print()
print("=== DEVOSIA TAXONOMIC SET ===")
print(f"Genomes classified as Devosia/Devosia_* : {len(devosia)}")

# ------------------------------------------------------------
# Non-Devosia
# ------------------------------------------------------------

non_dev = [
    r for r in integrated
    if not r["GTDB_Genus"].lower().startswith("devosia")
]

print()
print("=== NON-DEVOSIA ===")
print(f"Non-Devosia genomes : {len(non_dev)}")

# ------------------------------------------------------------
# ANI classification
# ------------------------------------------------------------

ani_only = [
    r for r in integrated
    if r["GTDB_Classification_Method"] == "ani_screen"
]

topology = [
    r for r in integrated
    if "topology" in r["GTDB_Note"].lower()
]

print()
print("=== CLASSIFICATION TYPE ===")
print(f"ANI-screen classifications : {len(ani_only)}")
print(f"Topology-based             : {len(topology)}")

# ------------------------------------------------------------
# Write Devosia-only table
# ------------------------------------------------------------

devosia_file = out_dir / "Devosia_GTDB_confirmed.tsv"

with open(devosia_file, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(devosia)

# ------------------------------------------------------------
# Write non-Devosia table
# ------------------------------------------------------------

nondev_file = out_dir / "Non_Devosia_GTDB.tsv"

with open(nondev_file, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(non_dev)

# ------------------------------------------------------------
# Write summary text
# ------------------------------------------------------------

summary_file = out_dir / "validation_summary.txt"

with open(summary_file, "w") as fh:

    fh.write("DEVOSIA GENOME VALIDATION SUMMARY\n")
    fh.write("=" * 70 + "\n\n")

    fh.write(f"CheckM2 genomes: {len(qc_rows)}\n")
    fh.write(f"GTDB-Tk genomes: {len(gtdb_rows)}\n")
    fh.write(f"Integrated genomes: {len(integrated)}\n\n")

    fh.write("CheckM2 QC STATUS\n")
    fh.write("-" * 40 + "\n")

    for k, v in qc_status.most_common():
        fh.write(f"{k}\t{v}\n")

    fh.write("\nGTDB GENUS\n")
    fh.write("-" * 40 + "\n")

    for k, v in genera.most_common():
        fh.write(f"{k}\t{v}\n")

    fh.write("\nGTDB FAMILY\n")
    fh.write("-" * 40 + "\n")

    for k, v in families.most_common():
        fh.write(f"{k}\t{v}\n")

    fh.write("\nCLASSIFICATION METHODS\n")
    fh.write("-" * 40 + "\n")

    for k, v in methods.most_common():
        fh.write(f"{k}\t{v}\n")

    fh.write("\n")
    fh.write(f"Species assigned: {species_assigned}\n")
    fh.write(f"Species unassigned: {len(integrated) - species_assigned}\n")
    fh.write(f"Devosia: {len(devosia)}\n")
    fh.write(f"Non-Devosia: {len(non_dev)}\n")
    fh.write(f"ANI-screen: {len(ani_only)}\n")
    fh.write(f"Topology-based: {len(topology)}\n")

print()
print("=" * 70)
print("FILES CREATED")
print("=" * 70)

print(integrated_file)
print(devosia_file)
print(nondev_file)
print(summary_file)

PY

echo
echo "============================================================"
echo "VALIDATION COMPLETE"
echo "============================================================"

echo
echo "=== Output files ==="
ls -lh "$OUT"

echo
echo "=== Integrated table ==="
column -ts $'\t' "$OUT/Devosia_GTDB_CheckM2_integrated.tsv" | head -10

echo
echo "=== Summary ==="
cat "$OUT/validation_summary.txt"
