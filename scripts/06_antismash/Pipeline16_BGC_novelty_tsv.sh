#!/usr/bin/env bash

# ============================================================
# Pipeline 16 — BGC novelty classification
#
# Purpose:
#   Classify antiSMASH BGC regions according to their maximum
#   ClusterBlast similarity to known MIBiG clusters.
#
# Categories:
#   novel   : max similarity < 30%
#   related : 30% <= max similarity < 70%
#   known   : max similarity >= 70%
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

## ---- Project paths ----

ANTISMASH_DIR="$PROJECT_ROOT/antismash_output"
OUTPUT="$PROJECT_ROOT/bgc_novelty_classification.tsv"

## ---- Configuration ----

NOVEL_MAX=30
KNOWN_MIN=70

## ---- Input validation ----

if [ ! -d "$ANTISMASH_DIR" ]; then
    echo "ERROR: antiSMASH output directory not found:"
    echo "  $ANTISMASH_DIR"
    exit 1
fi

JSON_COUNT=$(find "$ANTISMASH_DIR" -mindepth 2 -maxdepth 2 \
    -type f -name "*.json" | wc -l)

if [ "$JSON_COUNT" -eq 0 ]; then
    echo "ERROR: No antiSMASH JSON files found:"
    echo "  $ANTISMASH_DIR"
    exit 1
fi

## ---- Report configuration ----

echo "============================================================"
echo "Pipeline 16 — BGC novelty classification"
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
echo "Novel threshold:"
echo "  < ${NOVEL_MAX}%"

echo "Related range:"
echo "  ${NOVEL_MAX}% to < ${KNOWN_MIN}%"

echo "Known threshold:"
echo "  >= ${KNOWN_MIN}%"

echo
echo "Output:"
echo "  $OUTPUT"

## ---- Classification ----

python3 - "$ANTISMASH_DIR" "$OUTPUT" "$NOVEL_MAX" "$KNOWN_MIN" <<'PY'
#!/usr/bin/env python3

import sys
import json
import csv
from pathlib import Path
from collections import Counter


ANTISMASH_DIR = Path(sys.argv[1])
OUTPUT = Path(sys.argv[2])

NOVEL_MAX = float(sys.argv[3])
KNOWN_MIN = float(sys.argv[4])


rows = []


json_files = sorted(ANTISMASH_DIR.glob("*/*.json"))

print(f"Processing {len(json_files)} antiSMASH JSON files...")


for json_file in json_files:

    acc = json_file.stem

    with json_file.open() as handle:
        data = json.load(handle)

    for record in data.get("records", []):

        contig = record.get("id", "")

        clusterblast_module = (
            record
            .get("modules", {})
            .get("antismash.modules.clusterblast", {})
        )

        knowncluster = clusterblast_module.get("knowncluster", {})

        results = knowncluster.get("results", [])

        # Index ClusterBlast results by region number.
        by_region = {
            result["region_number"]: result
            for result in results
            if "region_number" in result
        }

        # Extract region numbers from region features.
        region_numbers = sorted({
            int(feature["qualifiers"]["region_number"][0])
            for feature in record.get("features", [])
            if (
                feature.get("type") == "region"
                and "region_number" in feature.get("qualifiers", {})
            )
        })

        for region_number in region_numbers:

            entry = by_region.get(region_number)

            if entry is None or not entry.get("ranking"):

                # No ClusterBlast hits:
                # classify as novel.
                max_similarity = 0
                n_hits = 0
                top_mibig = ""

            else:

                ranking = entry["ranking"]

                similarities = [
                    hit_info.get("similarity", 0)
                    for _, hit_info in ranking
                ]

                max_similarity = max(similarities)

                n_hits = len(ranking)

                top_index = similarities.index(max_similarity)

                top_hit = ranking[top_index][0]

                top_mibig = top_hit.get("accession", "")


            # ------------------------------------------------
            # Novelty classification
            # ------------------------------------------------

            if max_similarity < NOVEL_MAX:

                category = "novel"

            elif max_similarity < KNOWN_MIN:

                category = "related"

            else:

                category = "known"


            rows.append({
                "accession": acc,
                "contig": contig,
                "region_number": region_number,
                "max_similarity": max_similarity,
                "n_mibig_hits": n_hits,
                "top_mibig_hit": top_mibig,
                "novelty_category": category,
            })


# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

fieldnames = [
    "accession",
    "contig",
    "region_number",
    "max_similarity",
    "n_mibig_hits",
    "top_mibig_hit",
    "novelty_category",
]


OUTPUT.parent.mkdir(parents=True, exist_ok=True)


with OUTPUT.open("w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(rows)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print()
print("============================================================")
print("NOVELTY CLASSIFICATION SUMMARY")
print("============================================================")

print(f"Total regions classified: {len(rows)}")

category_counts = Counter(
    row["novelty_category"]
    for row in rows
)

total = len(rows)

for category in ["novel", "related", "known"]:

    n = category_counts.get(category, 0)

    percentage = (
        100 * n / total
        if total > 0
        else 0
    )

    print(
        f"{category}: "
        f"{n} "
        f"({percentage:.2f}%)"
    )

print()
print(f"Output written to:")
print(f"  {OUTPUT}")
PY

## ---- Output validation ----

if [ ! -s "$OUTPUT" ]; then
    echo
    echo "ERROR: Novelty classification output was not created or is empty."
    exit 1
fi

echo
echo "============================================================"
echo "Pipeline 16 completed successfully."
echo "============================================================"

echo
echo "Output:"
echo "  $OUTPUT"

echo
echo "Preview:"
head -10 "$OUTPUT"

echo
echo "Output line count:"
wc -l "$OUTPUT"

echo
