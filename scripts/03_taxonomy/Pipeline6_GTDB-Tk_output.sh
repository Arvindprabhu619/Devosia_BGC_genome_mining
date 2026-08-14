#!/usr/bin/env bash

# ============================================================
# Pipeline 6 — GTDB-Tk summary diagnostic
#
# Purpose:
#   Inspect the GTDB-Tk R232 bac120 classification output,
#   including the summary header, representative rows,
#   failed genomes, and available TSV outputs.
#
# Input:
#   gtdbtk_r232/classify/gtdbtk.bac120.summary.tsv
#
# Output:
#   Diagnostic information printed to stdout.
#
# Portability:
#   Repository root is inferred from this script's location.
# ============================================================

set -euo pipefail

## ---- Determine repository root ----

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

## ---- Input ----

FILE="$PROJECT_ROOT/gtdbtk_r232/classify/gtdbtk.bac120.summary.tsv"
FAILED_FILE="$PROJECT_ROOT/gtdbtk_r232/identify/gtdbtk.failed_genomes.tsv"

if [ ! -f "$FILE" ]; then
    echo "ERROR: GTDB-Tk summary file not found:"
    echo "  $FILE"
    echo
    echo "Run Pipeline4_GTDBTK_MAG.sh first."
    exit 1
fi

echo "============================================================"
echo "GTDB-Tk SUMMARY DIAGNOSTIC"
echo "============================================================"

echo
echo "=== Project root ==="
echo "$PROJECT_ROOT"

echo
echo "=== Summary file ==="
ls -lh "$FILE"

echo
echo "=== Number of lines ==="
wc -l "$FILE"

echo
echo "=== EXACT HEADER ==="
head -1 "$FILE" | cat -A

echo
echo "=== HEADER COLUMNS ==="

python3 - "$FILE" <<'PY'
import csv
import sys

f = sys.argv[1]

with open(f, newline="") as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)

print("Number of columns:", len(header))

for i, col in enumerate(header, 1):
    print(f"{i:2d}: {col!r}")
PY

echo
echo "=== FIRST 5 DATA ROWS ==="
head -6 "$FILE" | cut -c1-1000

echo
echo "=== IDENTIFY OUTPUT ==="

if [ -f "$FAILED_FILE" ]; then

    wc -l "$FAILED_FILE"

    echo
    cat "$FAILED_FILE"

else

    echo "No GTDB-Tk failed-genomes file found:"
    echo "  $FAILED_FILE"

fi

echo
echo "=== ALL GTDB-Tk OUTPUT SUMMARY FILES ==="

find "$PROJECT_ROOT/gtdbtk_r232" \
    -type f \
    -name "*.tsv" \
    -printf "%p\n" |
    sort

echo
echo "============================================================"
echo "GTDB-Tk diagnostic completed."
echo "============================================================"
