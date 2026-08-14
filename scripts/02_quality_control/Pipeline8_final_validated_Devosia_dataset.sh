
conda activate gtdbtk

export GTDBTK_DATA_PATH="$HOME/Microbial_minning_project/gtdbtk_db/release232"

BASE="$PROJECT_ROOT"
QC="$BASE/Devosia_CheckM2_QC.tsv"
GTDB="$BASE/gtdbtk_r232/classify/gtdbtk.bac120.summary.tsv"
GENOMES="$BASE/genomes"
OUT="$BASE/final_validation"

mkdir -p "$OUT"

echo "============================================================"
echo "FINAL DEVOSIA GENOME VALIDATION"
echo "CheckM2 + GTDB-Tk R232"
echo "============================================================"

date

python - "$QC" "$GTDB" "$GENOMES" "$OUT" <<'PY'

import csv
import sys
import shutil
from pathlib import Path
from collections import Counter

QC_FILE = Path(sys.argv[1])
GTDB_FILE = Path(sys.argv[2])
GENOMES_DIR = Path(sys.argv[3])
OUT_DIR = Path(sys.argv[4])

# ------------------------------------------------------------
# Read CheckM2
# ------------------------------------------------------------

with open(QC_FILE, newline="") as f:
    qc_rows = list(csv.DictReader(f, delimiter="\t"))

# ------------------------------------------------------------
# Read GTDB-Tk
# ------------------------------------------------------------

with open(GTDB_FILE, newline="") as f:
    gtdb_rows = list(csv.DictReader(f, delimiter="\t"))

qc = {r["Name"]: r for r in qc_rows}
gtdb = {r["user_genome"]: r for r in gtdb_rows}

print()
print("=" * 70)
print("INPUT VALIDATION")
print("=" * 70)

print(f"CheckM2 genomes : {len(qc)}")
print(f"GTDB-Tk genomes : {len(gtdb)}")

common = sorted(set(qc) & set(gtdb))

print(f"Matched genomes : {len(common)}")

if len(common) != len(qc) or len(common) != len(gtdb):
    print()
    print("WARNING: CheckM2 and GTDB-Tk genome sets differ.")

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def taxonomy_field(row, prefix):
    value = row["classification"]
    for x in value.split(";"):
        if x.startswith(prefix):
            return x.replace(prefix, "")
    return ""

def is_devosia(row):
    genus = taxonomy_field(row, "g__")
    return genus in {"Devosia", "Devosia_A"}

def is_devosia_family(row):
    family = taxonomy_field(row, "f__")
    return family == "Devosiaceae"

def species_name(row):
    return taxonomy_field(row, "s__")

# ------------------------------------------------------------
# CheckM2 retention
#
# Use the QC_status already assigned by the QC workflow.
# ------------------------------------------------------------

retained_statuses = {
    "RETAINED_STANDARD",
    "RETAINED_BORDERLINE"
}

# ------------------------------------------------------------
# Build integrated final table
# ------------------------------------------------------------

integrated = []

for genome in common:

    q = qc[genome]
    g = gtdb[genome]

    genus = taxonomy_field(g, "g__")
    family = taxonomy_field(g, "f__")
    species = species_name(g)

    retained_qc = q["QC_status"] in retained_statuses
    devosia_genus = is_devosia(g)
    devosia_family = is_devosia_family(g)

    if retained_qc and devosia_genus:
        final_status = "FINAL_RETAIN"
    elif not devosia_genus:
        final_status = "NON_DEVOSIA"
    elif not retained_qc:
        final_status = "CHECKM2_EXCLUDED"
    else:
        final_status = "REVIEW"

    integrated.append({
        "Genome": genome,
        "CheckM2_Completeness": q["Completeness"],
        "CheckM2_Contamination": q["Contamination"],
        "CheckM2_QC_status": q["QC_status"],

        "GTDB_Domain": taxonomy_field(g, "d__"),
        "GTDB_Phylum": taxonomy_field(g, "p__"),
        "GTDB_Class": taxonomy_field(g, "c__"),
        "GTDB_Order": taxonomy_field(g, "o__"),
        "GTDB_Family": family,
        "GTDB_Genus": genus,
        "GTDB_Species": species,

        "GTDB_Classification_Method":
            g["classification_method"],

        "GTDB_Note":
            g["note"],

        "Closest_Genome":
            g["closest_genome_reference"],

        "Closest_Genome_ANI":
            g["closest_genome_ani"],

        "Closest_Genome_AF":
            g["closest_genome_af"],

        "MSA_Percent":
            g["msa_percent"],

        "Translation_Table":
            g["translation_table"],

        "RED_Value":
            g["red_value"],

        "Warnings":
            g["warnings"],

        "Final_Status":
            final_status
    })

# ------------------------------------------------------------
# Write integrated table
# ------------------------------------------------------------

fields = list(integrated[0].keys())

integrated_file = OUT_DIR / "FINAL_Devosia_Validation.tsv"

with open(integrated_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(integrated)

# ------------------------------------------------------------
# Final retained genomes
# ------------------------------------------------------------

final_rows = [
    r for r in integrated
    if r["Final_Status"] == "FINAL_RETAIN"
]

final_file = OUT_DIR / "FINAL_Devosia_Genomes.tsv"

with open(final_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(final_rows)

# ------------------------------------------------------------
# Excluded genomes
# ------------------------------------------------------------

excluded_rows = [
    r for r in integrated
    if r["Final_Status"] == "CHECKM2_EXCLUDED"
]

excluded_file = OUT_DIR / "CheckM2_Excluded.tsv"

with open(excluded_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(excluded_rows)

# ------------------------------------------------------------
# Non-Devosia
# ------------------------------------------------------------

non_devosia_rows = [
    r for r in integrated
    if r["Final_Status"] == "NON_DEVOSIA"
]

non_devosia_file = OUT_DIR / "Non_Devosia_Genomes.tsv"

with open(non_devosia_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(non_devosia_rows)

# ------------------------------------------------------------
# Copy final FASTA files
# ------------------------------------------------------------

FINAL_FASTA = OUT_DIR / "final_genomes"

if FINAL_FASTA.exists():
    shutil.rmtree(FINAL_FASTA)

FINAL_FASTA.mkdir()

copied = 0
missing = []

for row in final_rows:

    genome = row["Genome"]

    source = GENOMES_DIR / f"{genome}.fna"

    if not source.exists():
        matches = list(GENOMES_DIR.glob(f"{genome}*.fna"))

        if matches:
            source = matches[0]
        else:
            missing.append(genome)
            continue

    shutil.copy2(source, FINAL_FASTA / source.name)
    copied += 1

# ------------------------------------------------------------
# Statistics
# ------------------------------------------------------------

status_counts = Counter(
    r["Final_Status"] for r in integrated
)

genus_counts = Counter(
    r["GTDB_Genus"]
    for r in final_rows
)

species_counts = Counter(
    r["GTDB_Species"]
    for r in final_rows
    if r["GTDB_Species"]
)

method_counts = Counter(
    r["GTDB_Classification_Method"]
    for r in final_rows
)

# ------------------------------------------------------------
# Summary report
# ------------------------------------------------------------

summary_file = OUT_DIR / "FINAL_validation_summary.txt"

with open(summary_file, "w") as f:

    f.write("=" * 70 + "\n")
    f.write("FINAL DEVOSIA GENOME VALIDATION\n")
    f.write("CheckM2 + GTDB-Tk R232\n")
    f.write("=" * 70 + "\n\n")

    f.write(f"CheckM2 genomes       : {len(qc)}\n")
    f.write(f"GTDB-Tk genomes       : {len(gtdb)}\n")
    f.write(f"Matched genomes       : {len(common)}\n\n")

    f.write("FINAL STATUS\n")
    f.write("-" * 40 + "\n")

    for k, v in status_counts.most_common():
        f.write(f"{k:30s} {v}\n")

    f.write("\n")
    f.write(f"FINAL RETAINED DEVOSIA : {len(final_rows)}\n")
    f.write(f"FINAL FASTA COPIED     : {copied}\n")
    f.write(f"FASTA MISSING          : {len(missing)}\n\n")

    f.write("GTDB-Tk GENUS DISTRIBUTION\n")
    f.write("-" * 40 + "\n")

    for k, v in genus_counts.most_common():
        f.write(f"{k:30s} {v}\n")

    f.write("\n")
    f.write("GTDB-Tk SPECIES DISTRIBUTION\n")
    f.write("-" * 40 + "\n")

    for k, v in species_counts.most_common():
        f.write(f"{k:50s} {v}\n")

    f.write("\n")
    f.write("GTDB-Tk CLASSIFICATION METHODS\n")
    f.write("-" * 40 + "\n")

    for k, v in method_counts.most_common():
        f.write(f"{k:50s} {v}\n")

    if missing:
        f.write("\n")
        f.write("MISSING FASTA FILES\n")
        f.write("-" * 40 + "\n")

        for x in missing:
            f.write(x + "\n")

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

print()
print("=" * 70)
print("FINAL VALIDATION RESULT")
print("=" * 70)

print()
print("=== FINAL STATUS ===")

for k, v in status_counts.most_common():
    print(f"{k:30s} {v}")

print()
print(f"FINAL RETAINED DEVOSIA : {len(final_rows)}")
print(f"FINAL FASTA COPIED     : {copied}")
print(f"FASTA MISSING          : {len(missing)}")

print()
print("=== FINAL GENUS DISTRIBUTION ===")

for k, v in genus_counts.most_common():
    print(f"{k:30s} {v}")

print()
print("=== FINAL SPECIES COUNT ===")

print(f"Unique species assigned : {len(species_counts)}")

print()
print("=== CLASSIFICATION METHODS ===")

for k, v in method_counts.most_common():
    print(f"{k:50s} {v}")

print()
print("============================================================")
print("OUTPUT FILES")
print("============================================================")

print(integrated_file)
print(final_file)
print(excluded_file)
print(non_devosia_file)
print(summary_file)
print(FINAL_FASTA)

PY

echo
echo "============================================================"
echo "FINAL VALIDATION DIRECTORY"
echo "============================================================"

find "$OUT" -maxdepth 2 -type f \
    -printf "%p\t%k KB\n" | sort

echo
echo "============================================================"
echo "FINAL GENOME COUNT"
echo "============================================================"

find "$OUT/final_genomes" \
    -maxdepth 1 \
    -type f \
    -name "*.fna" | wc -l

echo
echo "============================================================"
echo "DONE"
echo "============================================================"
