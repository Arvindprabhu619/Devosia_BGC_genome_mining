#!/usr/bin/env bash

# ============================================================
# Pipeline 4 — GTDB-Tk R232 classification
#
# Purpose:
#   Classify the curated Devosia genome dataset using GTDB-Tk.
#
# External dependency:
#   GTDB-Tk R232 reference database.
#
# Before running:
#   export GTDBTK_DATA_PATH=/path/to/gtdbtk/release232
#
# Inputs:
#   genomes/
#
# Output:
#   gtdbtk_r232/
# ============================================================

set -euo pipefail

## ---- Determine repository root ----

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

## ---- Validate GTDB-Tk database configuration ----

if [ -z "${GTDBTK_DATA_PATH:-}" ]; then
    echo "ERROR: GTDBTK_DATA_PATH is not set."
    echo
    echo "Set it to your local GTDB-Tk R232 database, for example:"
    echo
    echo "  export GTDBTK_DATA_PATH=/path/to/gtdbtk_db/release232"
    echo
    exit 1
fi

if [ ! -d "$GTDBTK_DATA_PATH" ]; then
    echo "ERROR: GTDB-Tk database directory does not exist:"
    echo "  $GTDBTK_DATA_PATH"
    exit 1
fi

## ---- Activate GTDB-Tk environment ----

if ! command -v gtdbtk >/dev/null 2>&1; then

    echo "ERROR: gtdbtk command not found."
    echo
    echo "Activate the appropriate GTDB-Tk environment first."
    echo "Example:"
    echo
    echo "  conda activate gtdbtk"
    echo

    exit 1

fi

## ---- Input/output ----

GENOMES="$PROJECT_ROOT/genomes"
GTDB_OUT="$PROJECT_ROOT/gtdbtk_r232"

if [ ! -d "$GENOMES" ]; then
    echo "ERROR: Genome directory not found:"
    echo "  $GENOMES"
    exit 1
fi

mkdir -p "$GTDB_OUT"

echo "============================================================"
echo "GTDB-Tk R232 CLASSIFICATION"
echo "============================================================"

echo "Project root:"
echo "  $PROJECT_ROOT"

echo "Genome directory:"
echo "  $GENOMES"

echo "GTDB-Tk database:"
echo "  $GTDBTK_DATA_PATH"

echo "Output directory:"
echo "  $GTDB_OUT"

echo

date

gtdbtk classify_wf \
    --genome_dir "$GENOMES" \
    --out_dir "$GTDB_OUT" \
    --cpus 32 \
    --force

STATUS=$?

echo
echo "============================================================"
echo "GTDB-Tk EXIT STATUS: $STATUS"
echo "============================================================"

date

if [ "$STATUS" -eq 0 ]; then

    echo "GTDB-Tk classification completed successfully."

else

    echo "GTDB-Tk classification FAILED."

    exit "$STATUS"

fi
