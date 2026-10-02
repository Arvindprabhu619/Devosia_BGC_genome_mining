#!/usr/bin/env bash
# ============================================================
# scripts/06_run_antismash.sh
# Purpose: Predict BGCs in each genome with antiSMASH 7.1.0
# Inputs:  results/04_gtdbtk/final_genome_list.txt
# Outputs: results/06_antismash/raw_output/<acc>/<acc>.gbk
# Runtime: 12-24 hours
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
mkdir -p "${RES_ANTISMASH}/logs"

# Prepare inputs
echo "[06.1] Preparing input files..."
n=0
while read -r acc; do
    [[ -z "$acc" ]] && continue
    fasta=$(awk -F'\t' -v a="$acc" '$1==a {print $2}' "${DATA_META}/genome_manifest.tsv")
    gff=$(awk -F'\t' -v a="$acc" '$1==a {print $3}' "${DATA_META}/genome_manifest.tsv")
    if [[ -z "$fasta" || ! -f "$fasta" ]]; then
        echo "  WARNING: no FASTA for $acc"
        continue
    fi
    ln -sf "$fasta" "${RES_ANTISMASH}/input/${acc}.fna"
    if [[ -n "$gff" && -f "$gff" ]]; then
        ln -sf "$gff" "${RES_ANTISMASH}/input/${acc}.gff"
    fi
    n=$((n+1))
done < "$FINAL_LIST"
echo "  Prepared $n genomes"

# Run function — log OUTSIDE the output dir
run_one() {
    local acc="$1"
    local fasta="${RES_ANTISMASH}/input/${acc}.fna"
    local outdir="${RES_ANTISMASH}/raw_output/${acc}"
    local logdir="${RES_ANTISMASH}/logs"

    # Completely remove anything from prior attempts
    rm -rf "$outdir"
    mkdir -p "$outdir"

    # Run antiSMASH — redirect log to SEPARATE logdir
    "$ANTI_SMASH" \
        --output-dir "$outdir" \
        --cpus 4 \
        --cb-general \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --rre \
        --genefinding-tool prodigal \
        --minlength 1000 \
        "$fasta" > "${logdir}/${acc}.log" 2>&1

    if [[ -f "${outdir}/${acc}.gbk" ]]; then
        echo "DONE: $acc"
    else
        echo "FAILED: $acc"
    fi
}

export -f run_one
export ANTI_SMASH RES_ANTISMASH

echo "[06.2] Running antiSMASH in parallel (12-24 hours)..."
cat "$FINAL_LIST" | xargs -P 4 -I {} bash -c 'run_one "$@"' _ {}

echo "[06.3] Summary:"
echo "  Successful: $(find ${RES_ANTISMASH}/raw_output -name "*.gbk" 2>/dev/null | wc -l)"
echo "  Failed:     $(find ${RES_ANTISMASH}/raw_output -type d 2>/dev/null | tail -n +2 | while read d; do [[ ! -f "$d"/*.gbk ]] && echo X; done | wc -l)"

echo "=== [06] Done: $(date) ==="
