#!/usr/bin/env bash

# ============================================================
# Pipeline 12a — GTDB-Tk phylogenomic tree generation
#
# Purpose:
#   Generate the final GTDB-Tk BAC120 phylogeny for the
#   curated Devosia genome dataset.
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

## ---- Project directories ----

PHYLO_DIR="$PROJECT_ROOT/phylogeny_110"

GENOME_DIR="$PHYLO_DIR/input_110_rooted"
TAXONOMY_FILE="$PHYLO_DIR/custom_taxonomy_115.tsv"
OUT_DIR="$PHYLO_DIR/gtdbtk_final_rooted_v2"

## ---- Input checks ----

if [ ! -d "$GENOME_DIR" ]; then
    echo "ERROR: Genome directory not found:"
    echo "  $GENOME_DIR"
    exit 1
fi

if [ ! -f "$TAXONOMY_FILE" ]; then
    echo "ERROR: Custom taxonomy file not found:"
    echo "  $TAXONOMY_FILE"
    exit 1
fi

## ---- Report configuration ----

echo "============================================================"
echo "Pipeline 12a — GTDB-Tk phylogenomic tree"
echo "============================================================"

echo
echo "Project root:"
echo "  $PROJECT_ROOT"

echo
echo "Genome directory:"
echo "  $GENOME_DIR"

echo
echo "Custom taxonomy:"
echo "  $TAXONOMY_FILE"

echo
echo "Output directory:"
echo "  $OUT_DIR"

echo
echo "Genome count:"
find "$GENOME_DIR" -maxdepth 1 -type f \
    \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) |
    wc -l

## ---- GTDB-Tk ----

gtdbtk de_novo_wf \
    --genome_dir "$GENOME_DIR" \
    --bacteria \
    --outgroup_taxon g__Paradevosia \
    --skip_gtdb_refs \
    --custom_taxonomy_file "$TAXONOMY_FILE" \
    --out_dir "$OUT_DIR" \
    --cpus 16

echo
echo "============================================================"
echo "Pipeline 12a completed."
echo "============================================================"

echo
echo "GTDB-Tk output:"
echo "  $OUT_DIR"
