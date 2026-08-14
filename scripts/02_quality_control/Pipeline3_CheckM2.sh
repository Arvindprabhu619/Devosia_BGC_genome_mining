#!/bin/bash

## ---- Determine repository root ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"


# ============================================================
# Devosia Genome Quality Assessment using CheckM2
# ============================================================
#
# Dataset:
#   137 curated Devosia genome assemblies
#
# Software:
#   CheckM2 1.1.0
#   Python 3.12.13
#   Conda 26.1.1
#
# CheckM2 database:
#   UniRef100 KEGG Orthology database
#
# QC criteria:
#
#   STANDARD:
#       Completeness >= 95%
#       Contamination <= 5.0%
#
#   BORDERLINE:
#       Completeness >= 95%
#       Contamination > 5.0% and <= 5.5%
#
#   EXCLUDED:
#       Completeness < 95%
#       OR contamination > 5.5%
#
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

BASE_DIR="$PROJECT_ROOT"

GENOMES_DIR="$BASE_DIR/genomes"

CHECKM2_OUT="$BASE_DIR/checkm2_out"

QC_TABLE="$BASE_DIR/Devosia_CheckM2_QC.tsv"

CHECKM2_DB="$HOME/databases/CheckM2_database/uniref100.KO.1.dmnd"

EXPECTED_GENOMES=137

# ------------------------------------------------------------
# CPU configuration
# ------------------------------------------------------------
#
# If running inside SLURM:
# use allocated CPUs.
#
# Otherwise:
# use available CPUs.
#

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
    THREADS="$SLURM_CPUS_PER_TASK"
else
    THREADS="$(nproc)"
fi

# ------------------------------------------------------------
# 2. Header
# ------------------------------------------------------------

echo
echo "============================================================"
echo "        Devosia Genome Quality Assessment - CheckM2"
echo "============================================================"
echo
echo "Base directory : $BASE_DIR"
echo "Genome directory: $GENOMES_DIR"
echo "Output directory: $CHECKM2_OUT"
echo "QC table        : $QC_TABLE"
echo "Database        : $CHECKM2_DB"
echo "Threads         : $THREADS"
echo

# ------------------------------------------------------------
# 3. Directory checks
# ------------------------------------------------------------

if [ ! -d "$BASE_DIR" ]; then
    echo "ERROR: Base directory does not exist:"
    echo "$BASE_DIR"
    exit 1
fi

if [ ! -d "$GENOMES_DIR" ]; then
    echo "ERROR: Genome directory does not exist:"
    echo "$GENOMES_DIR"
    exit 1
fi

# ------------------------------------------------------------
# 4. Check CheckM2
# ------------------------------------------------------------

if ! command -v checkm2 >/dev/null 2>&1; then
    echo "ERROR: CheckM2 command not found."
    echo "Run: conda activate checkm2"
    exit 1
fi

echo "=== SOFTWARE VERSIONS ==="

echo "CheckM2:"
checkm2 --version

echo

echo "Python:"
python --version

echo

echo "Conda:"
conda --version

echo

# ------------------------------------------------------------
# 5. Verify genome count
# ------------------------------------------------------------

echo "=== GENOME COUNT CHECK ==="

N_GENOMES=$(find "$GENOMES_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*.fna" \
    | wc -l)

echo "Expected genomes : $EXPECTED_GENOMES"
echo "Found genomes    : $N_GENOMES"

if [ "$N_GENOMES" -ne "$EXPECTED_GENOMES" ]; then

    echo
    echo "ERROR: Expected $EXPECTED_GENOMES genomes but found $N_GENOMES."
    echo
    echo "Genome files:"
    
    find "$GENOMES_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.fna" \
        -printf "%f\n" \
        | sort

    exit 1
fi

echo
echo "Genome count verified: $EXPECTED_GENOMES"
echo

# ------------------------------------------------------------
# 6. Check for empty genomes
# ------------------------------------------------------------

echo "=== GENOME FILE INTEGRITY ==="

EMPTY_FILES=$(find "$GENOMES_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*.fna" \
    -size 0)

if [ -n "$EMPTY_FILES" ]; then

    echo
    echo "ERROR: Empty genome files detected:"
    echo "$EMPTY_FILES"

    exit 1
fi

echo "No empty genome files detected."
echo

# ------------------------------------------------------------
# 7. Verify CheckM2 database
# ------------------------------------------------------------

echo "=== DATABASE CHECK ==="

if [ ! -f "$CHECKM2_DB" ]; then

    echo
    echo "ERROR: CheckM2 database not found:"
    echo "$CHECKM2_DB"
    echo

    echo "Current CheckM2 database configuration:"
    checkm2 database --current

    exit 1
fi

echo "Database found:"
echo "$CHECKM2_DB"

echo
echo "Database size:"
ls -lh "$CHECKM2_DB"

echo
echo "CheckM2 database registration:"
checkm2 database --current

echo

# ------------------------------------------------------------
# 8. Check existing output
# ------------------------------------------------------------

echo "=== EXISTING OUTPUT CHECK ==="

if [ -d "$CHECKM2_OUT" ]; then

    echo
    echo "WARNING: Output directory already exists:"
    echo "$CHECKM2_OUT"
    echo

    if [ -f "$CHECKM2_OUT/quality_report.tsv" ]; then

        N_EXISTING=$(tail -n +2 \
            "$CHECKM2_OUT/quality_report.tsv" \
            | wc -l)

        echo "Existing quality report:"
        echo "$N_EXISTING genome records."

        if [ "$N_EXISTING" -eq "$EXPECTED_GENOMES" ]; then

            echo
            echo "A complete CheckM2 result already exists."
            echo
            echo "No new calculation will be performed."
            echo
            echo "Existing report:"
            echo "$CHECKM2_OUT/quality_report.tsv"

            exit 0

        fi

    fi

    echo
    echo "Existing output is incomplete or invalid."
    echo
    echo "If you want a clean rerun, remove it with:"
    echo
    echo "rm -rf \"$CHECKM2_OUT\""
    echo
    echo "Then execute this script again."

    exit 1

fi

# ------------------------------------------------------------
# 9. Run CheckM2
# ------------------------------------------------------------

echo
echo "============================================================"
echo "                    RUNNING CHECKM2"
echo "============================================================"
echo

echo "Input genomes : $GENOMES_DIR"
echo "Output        : $CHECKM2_OUT"
echo "Database      : $CHECKM2_DB"
echo "Threads       : $THREADS"
echo

echo "Start time:"
date

echo

checkm2 predict \
    --input "$GENOMES_DIR" \
    --output-directory "$CHECKM2_OUT" \
    --extension fna \
    --threads "$THREADS" \
    --database_path "$CHECKM2_DB"

echo

echo "End time:"
date

echo

echo "============================================================"
echo "                 CHECKM2 RUN COMPLETED"
echo "============================================================"
echo

# ------------------------------------------------------------
# 10. Validate CheckM2 output
# ------------------------------------------------------------

QUALITY_REPORT="$CHECKM2_OUT/quality_report.tsv"

if [ ! -f "$QUALITY_REPORT" ]; then

    echo "ERROR: quality_report.tsv was not generated."

    echo
    echo "Check output directory:"
    echo "$CHECKM2_OUT"

    exit 1
fi

N_EVALUATED=$(tail -n +2 \
    "$QUALITY_REPORT" \
    | wc -l)

echo "Expected genomes : $EXPECTED_GENOMES"
echo "Evaluated        : $N_EVALUATED"

if [ "$N_EVALUATED" -ne "$EXPECTED_GENOMES" ]; then

    echo
    echo "ERROR: CheckM2 did not evaluate all 137 genomes."

    exit 1
fi

echo
echo "SUCCESS: All 137 genomes were evaluated."
echo

# ------------------------------------------------------------
# 11. Display CheckM2 report header
# ------------------------------------------------------------

echo "=== CHECKM2 QUALITY REPORT ==="

head -n 3 "$QUALITY_REPORT"

echo

# ------------------------------------------------------------
# 12. Generate independent QC classification
# ------------------------------------------------------------

echo "=== GENERATING QC CLASSIFICATION ==="

python << PYTHON

import csv
from pathlib import Path

quality_report = Path("$QUALITY_REPORT")
qc_output = Path("$QC_TABLE")

with quality_report.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    fieldnames = reader.fieldnames

if not fieldnames:
    raise RuntimeError("Could not read CheckM2 quality_report.tsv")

print("Columns detected:")
for column in fieldnames:
    print("  ", column)

# --------------------------------------------------------
# Identify completeness and contamination columns
# --------------------------------------------------------

def find_column(fieldnames, keyword):

    for column in fieldnames:

        if keyword.lower() in column.lower():
            return column

    return None


completeness_col = find_column(
    fieldnames,
    "completeness"
)

contamination_col = find_column(
    fieldnames,
    "contamination"
)

if completeness_col is None:
    raise RuntimeError(
        "Completeness column not found."
    )

if contamination_col is None:
    raise RuntimeError(
        "Contamination column not found."
    )

print()
print("Completeness column :", completeness_col)
print("Contamination column:", contamination_col)

# --------------------------------------------------------
# Classification
# --------------------------------------------------------

for row in rows:

    completeness = float(
        row[completeness_col]
    )

    contamination = float(
        row[contamination_col]
    )

    if (
        completeness >= 95.0
        and contamination <= 5.0
    ):

        status = "RETAINED_STANDARD"

    elif (
        completeness >= 95.0
        and contamination > 5.0
        and contamination <= 5.5
    ):

        status = "RETAINED_BORDERLINE"

    else:

        status = "EXCLUDED"

    row["QC_status"] = status


# --------------------------------------------------------
# Write QC table
# --------------------------------------------------------

output_fields = fieldnames + ["QC_status"]

with qc_output.open(
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=output_fields,
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(rows)


# --------------------------------------------------------
# Summary
# --------------------------------------------------------

standard = sum(
    row["QC_status"] == "RETAINED_STANDARD"
    for row in rows
)

borderline = sum(
    row["QC_status"] == "RETAINED_BORDERLINE"
    for row in rows
)

excluded = sum(
    row["QC_status"] == "EXCLUDED"
    for row in rows
)

retained = standard + borderline

print()
print("================================================")
print("                 QC SUMMARY")
print("================================================")

print(
    f"Total genomes       : {len(rows)}"
)

print(
    f"Standard retained   : {standard}"
)

print(
    f"Borderline retained : {borderline}"
)

print(
    f"Excluded            : {excluded}"
)

print(
    f"Total retained      : {retained}"
)

print("================================================")

PYTHON

echo

# ------------------------------------------------------------
# 13. Display borderline genomes
# ------------------------------------------------------------

echo "=== BORDERLINE GENOMES ==="

awk -F'\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if(tolower($i) ~ /completeness/)
            comp=i
        if(tolower($i) ~ /contamination/)
            cont=i
        if($i=="QC_status")
            status=i
    }
    print
    next
}
$status=="RETAINED_BORDERLINE" {
    print
}
' "$QC_TABLE"

echo

# ------------------------------------------------------------
# 14. Display excluded genomes
# ------------------------------------------------------------

echo "=== EXCLUDED GENOMES ==="

awk -F'\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="QC_status")
            status=i
    }
    print
    next
}
$status=="EXCLUDED" {
    print
}
' "$QC_TABLE"

echo

# ------------------------------------------------------------
# 15. Final output
# ------------------------------------------------------------

echo "============================================================"
echo "                     FINAL OUTPUT"
echo "============================================================"
echo

echo "Raw CheckM2 report:"
echo "  $QUALITY_REPORT"

echo

echo "QC classification table:"
echo "  $QC_TABLE"

echo

echo "CheckM2 output directory:"
echo "  $CHECKM2_OUT"

echo

echo "Completed:"
date

echo
echo "============================================================"
echo "                  PIPELINE FINISHED"
echo "============================================================"
