#!/usr/bin/env bash
# ============================================================
# scripts/03_build_manifest.sh
# ============================================================
# Purpose: Build a manifest mapping each accession to its FASTA/GFF3/FAA files.
#
# Inputs:   data/metadata/curated_accessions.txt
#           data/raw_genomes/ncbi_dataset/data/
# Outputs:
#   data/metadata/genome_manifest.tsv
#   data/metadata/genomes_with_protein.txt
#   data/metadata/genomes_without_protein.txt
#
# Runtime:  2 min
# ============================================================
set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/03_manifest_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [03] Build manifest: $(date) ==="

RAW="${DATA_RAW}/ncbi_dataset/data"
META="${DATA_META}"

if [[ ! -d "$RAW" ]]; then
    echo "ERROR: Raw genomes not found at $RAW"
    exit 1
fi

MANIFEST="${META}/genome_manifest.tsv"
echo -e "accession\tfasta\tgff\tfaa\thas_protein" > "$MANIFEST"

WITH=0
WITHOUT=0
> "${META}/genomes_with_protein.txt"
> "${META}/genomes_without_protein.txt"

for dir in "$RAW"/*/; do
    [[ ! -d "$dir" ]] && continue
    acc=$(basename "$dir")

    # Find genomic FASTA (not cds_from_genomic)
    fasta=$(find "$dir" -maxdepth 1 -name "*.fna" -not -name "cds_from_genomic.fna" 2>/dev/null | head -1)
    gff=$(find "$dir" -maxdepth 1 -name "*.gff" 2>/dev/null | head -1)
    faa=$(find "$dir" -maxdepth 1 -name "*.faa" 2>/dev/null | head -1)

    if [[ -z "$fasta" ]]; then
        continue
    fi

    has_protein="no"
    if [[ -n "$gff" && -n "$faa" ]]; then
        has_protein="yes"
        WITH=$((WITH + 1))
        echo "$acc" >> "${META}/genomes_with_protein.txt"
    else
        WITHOUT=$((WITHOUT + 1))
        echo "$acc" >> "${META}/genomes_without_protein.txt"
    fi

    echo -e "${acc}\t${fasta}\t${gff}\t${faa}\t${has_protein}" >> "$MANIFEST"
done

total=$(tail -n +2 "$MANIFEST" | wc -l)
echo ""
echo "=== Manifest summary ==="
echo "Total genomes:      $total"
echo "With protein/GFF3:  $WITH"
echo "Without protein:    $WITHOUT"

echo "=== [03] Done: $(date) ==="
