#!/usr/bin/env bash
# ============================================================
# scripts/06_run_antismash.sh
#
# Purpose: Predict biosynthetic gene clusters (BGCs) in each genome.
# Inputs:  results/04_gtdbtk/final_genome_list.txt (Devosia + Devosia_A)
#          Genome FASTA + GFF3 from manifest
# Outputs: results/06_antismash/raw_output/<accession>/  (per-genome)
#          results/06_antismash/summary.txt
#
# Tool: antiSMASH 7.1.0
# Runtime: 12-24 hours (with 16 parallel genomes)
# ============================================================

set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/06_antismash_$(date +%Y%m%d_%H%M).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [06] antiSMASH: $(date) ==="

FINAL_LIST="${RES_GTDBTK}/final_genome_list.txt"
if [[ ! -f "$FINAL_LIST" ]]; then
    echo "ERROR: final_genome_list.txt not found. Run GTDB-Tk first."
    exit 1
fi

mkdir -p "${RES_ANTISMASH}/raw_output"
mkdir -p "${RES_ANTISMASH}/input"

# Prepare input: copy FASTA + GFF3 for each final genome
echo "[06.1] Preparing input files..."
n=0
while read -r acc; do
    [[ -z "$acc" ]] && continue
    # Look up paths in manifest
    fasta=$(awk -F'\t' -v a="$acc" '$1==a {print $2}' "${DATA_META}/genome_manifest.tsv")
    gff=$(awk -F'\t' -v a="$acc" '$1==a {print $3}' "${DATA_META}/genome_manifest.tsv")

    if [[ -z "$fasta" || ! -f "$fasta" ]]; then
        echo "  WARNING: no FASTA for $acc, skipping"
        continue
    fi

    ln -sf "$fasta" "${RES_ANTISMASH}/input/${acc}.fna"
    if [[ -n "$gff" && -f "$gff" ]]; then
        ln -sf "$gff" "${RES_ANTISMASH}/input/${acc}.gff"
    fi
    n=$((n+1))
done < "$FINAL_LIST"
echo "  Prepared $n genomes"

# Run antiSMASH in parallel
echo "[06.2] Running antiSMASH in parallel (12-24 hours)..."

run_one() {
    local acc="$1"
    local fasta="${RES_ANTISMASH}/input/${acc}.fna"
    local outdir="${RES_ANTISMASH}/raw_output/${acc}"

    [[ -d "$outdir" ]] && return 0

    mkdir -p "$outdir"

    $ANTI_SMASH \
        --output-dir "$outdir" \
        --cpus 4 \
        --cb-general \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --rre \
        --smcogs \
        --genefinding-tool prodigal \
        --minlength 1000 \
        "$fasta" > "${outdir}/antismash.stdout.log" 2>&1

    if [[ -f "${outdir}/${acc}.gbk" ]]; then
        echo "DONE: $acc"
    else
        echo "FAILED: $acc"
    fi
}

export -f run_one
export ANTI_SMASH RES_ANTISMASH

# Run 4 genomes in parallel (each uses 4 CPUs = 16 total)
cat "$FINAL_LIST" | xargs -P 4 -I {} bash -c 'run_one "$@"' _ {}

echo "[06.3] antiSMASH complete. Summary:"
echo "  Successful: $(find ${RES_ANTISMASH}/raw_output -name "*.gbk" | wc -l)"
echo "  Failed:     $(find ${RES_ANTISMASH}/raw_output -name "antismash.stdout.log" -exec grep -L "antiSMASH" {} \; | wc -l)"

echo "=== [06] Done: $(date) ==="
