#!/usr/bin/env bash

# ============================================================
# Pipeline 14 — Extract antiSMASH BGC regions and classes
#
# Purpose:
#   Parse antiSMASH JSON files and generate:
#
#     bgc_regions_detailed.tsv
#
#   containing one row per predicted BGC region.
#
#   Hybrid/multi-class regions retain all antiSMASH product
#   annotations separated by semicolons.
#
# Portability:
#   Repository root is inferred from this script's location.
#   No hard-coded HOME or project-specific paths.
# ============================================================

set -euo pipefail

## ---- Determine repository root ----

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

## ---- Input/output ----

ANTISMASH_DIR="$PROJECT_ROOT/antismash_output"
OUTPUT="$PROJECT_ROOT/bgc_regions_detailed.tsv"

## ---- Input validation ----

if [ ! -d "$ANTISMASH_DIR" ]; then
    echo "ERROR: antiSMASH output directory not found:"
    echo "  $ANTISMASH_DIR"
    exit 1
fi

JSON_COUNT=$(find "$ANTISMASH_DIR" -mindepth 2 -maxdepth 2 \
    -type f -name "*.json" | wc -l)

if [ "$JSON_COUNT" -eq 0 ]; then
    echo "ERROR: No antiSMASH JSON files found in:"
    echo "  $ANTISMASH_DIR"
    exit 1
fi

## ---- Report configuration ----

echo "============================================================"
echo "Pipeline 14 — antiSMASH BGC extraction"
echo "============================================================"

echo
echo "Project root:"
echo "  $PROJECT_ROOT"

echo
echo "antiSMASH directory:"
echo "  $ANTISMASH_DIR"

echo
echo "JSON files detected:"
echo "  $JSON_COUNT"

echo
echo "Output:"
echo "  $OUTPUT"

## ---- Parse antiSMASH JSON ----

python3 - "$ANTISMASH_DIR" "$OUTPUT" <<'PY'
import sys
import json
import csv
from pathlib import Path
from collections import Counter

ANTISMASH_DIR = Path(sys.argv[1])
OUTPUT = Path(sys.argv[2])

rows = []

json_files = sorted(
    ANTISMASH_DIR.glob("*/*.json")
)

for f in json_files:

    acc = f.stem

    with open(f) as handle:
        data = json.load(handle)

    for record in data.get("records", []):

        contig = record.get("id", "")

        for area in record.get("areas", []):

            products = area.get("products", [])

            rows.append({
                "accession": acc,
                "contig": contig,
                "products": ";".join(products) if products else "unknown"
            })


fields = [
    "accession",
    "contig",
    "products"
]

with open(OUTPUT, "w", newline="") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)


print()
print("============================================================")
print("BGC EXTRACTION SUMMARY")
print("============================================================")

print(f"antiSMASH JSON files : {len(json_files)}")
print(f"BGC regions extracted: {len(rows)}")


# ------------------------------------------------------------
# Class-level tally
#
# A hybrid region can contain multiple product tags.
# Therefore class counts represent class-tag instances and
# can exceed the number of distinct BGC regions.
# ------------------------------------------------------------

class_counts = Counter()

for row in rows:

    for product in row["products"].split(";"):

        class_counts[product] += 1


print()
print("--- BGC class frequency ---")

for cls, count in class_counts.most_common():

    print(f"{cls}\t{count}")


print()
print("Output:")
print(OUTPUT)

PY

echo
echo "============================================================"
echo "Pipeline 14 completed."
echo "============================================================"

echo
echo "Output file:"
ls -lh "$OUTPUT"

echo
echo "First 10 records:"
head -10 "$OUTPUT"
