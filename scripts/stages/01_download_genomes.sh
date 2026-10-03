#!/usr/bin/env bash
# ============================================================
# scripts/01_download_genomes.sh
# ============================================================
# Purpose: Query NCBI for Devosia assemblies and download FASTA/GFF3/FAA.
#
# Inputs:   NCBI public database
# Outputs:
#   data/metadata/devosia_refseq.jsonl
#   data/metadata/devosia_genbank.jsonl
#   data/metadata/all_accessions.txt
#   data/raw_genomes/ncbi_dataset/
#
# Tool:     NCBI Datasets CLI 18.38.0
# Runtime:  15-30 min (depends on NCBI server speed)
# ============================================================
set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/01_download_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [01] Download genomes: $(date) ==="

cd "$DATA_META"

# ---- Query NCBI ----
if [[ ! -s "devosia_refseq.jsonl" ]]; then
    echo "[01.1] Querying NCBI for Devosia RefSeq assemblies..."
    $DATASETS summary genome taxon "Devosia" --assembly-source RefSeq \
        --as-json-lines > devosia_refseq.jsonl
    echo "  RefSeq: $(wc -l < devosia_refseq.jsonl) records"
else
    echo "[01.1] RefSeq JSONL already exists, skipping query"
fi

if [[ ! -s "devosia_genbank.jsonl" ]]; then
    echo "[01.2] Querying NCBI for Devosia GenBank assemblies..."
    $DATASETS summary genome taxon "Devosia" --assembly-source GenBank \
        --as-json-lines > devosia_genbank.jsonl
    echo "  GenBank: $(wc -l < devosia_genbank.jsonl) records"
else
    echo "[01.2] GenBank JSONL already exists, skipping query"
fi

# ---- Download ----
cd "$DATA_RAW"

if [[ ! -d "ncbi_dataset/data" ]] || [[ $(find ncbi_dataset/data -name "*.fna" 2>/dev/null | wc -l) -lt 100 ]]; then
    echo "[01.3] Downloading all Devosia genomes (10-30 min)..."
    $DATASETS download genome taxon "Devosia" \
        --include genome,gff3,protein \
        --filename "${DATA_META}/ncbi_dataset.zip" \
        --no-progressbar

    echo "[01.4] Extracting archive..."
    unzip -q -o "${DATA_META}/ncbi_dataset.zip" -d "${DATA_RAW}/ncbi_dataset"
else
    echo "[01.3] Genomes already downloaded, skipping"
fi

# ---- Summary ----
echo "[01.5] Summary:"
echo "  Genomes:  $(find ncbi_dataset/data -name "*.fna" 2>/dev/null | wc -l)"
echo "  GFF3:     $(find ncbi_dataset/data -name "*.gff" 2>/dev/null | wc -l)"
echo "  Protein:  $(find ncbi_dataset/data -name "*.faa" 2>/dev/null | wc -l)"

echo "=== [01] Done: $(date) ==="
