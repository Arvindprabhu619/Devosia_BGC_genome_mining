#!/usr/bin/env bash
# ============================================================
# scripts/09_run_bigscape.sh
#
# Purpose: Cluster BGCs into Gene Cluster Families (GCFs) using BiG-SCAPE.
# Inputs:  results/06_antismash/raw_output/<accession>/<accession>.gbk
# Outputs: results/07_bigscape/raw_output/  (BiG-SCAPE networks)
#          results/07_bigscape/parsed/gcf_assignments_c{0.3,0.5,0.7}.tsv
#
# Tool: BiG-SCAPE 2.0.3
# Runtime: 2-4 hours
# ============================================================

set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/09_bigscape_$(date +%Y%m%d_%H%M).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [09] BiG-SCAPE: $(date) ==="
echo "Pfam:      $PFAM_HMM"
echo "MIBiG dir: $MIBIG_JSON_DIR"

# Verify inputs
if [[ ! -f "$PFAM_HMM" ]]; then
    echo "ERROR: Pfam-A.hmm not found at $PFAM_HMM"
    exit 1
fi

mkdir -p "${RES_BIGSCAPE}/input"
mkdir -p "${RES_BIGSCAPE}/raw_output"
mkdir -p "${RES_BIGSCAPE}/parsed"

# Collect all antiSMASH GenBank files (one per genome)
echo "[09.1] Collecting antiSMASH GenBank files..."
n=0
for gbk in "${RES_ANTISMASH}"/raw_output/*/*.gbk; do
    [[ -f "$gbk" ]] || continue
    acc=$(basename "$(dirname "$gbk")")
    ln -sf "$gbk" "${RES_BIGSCAPE}/input/${acc}.gbk"
    n=$((n+1))
done
echo "  Linked $n GenBank files"

# Run BiG-SCAPE
echo "[09.2] Running BiG-SCAPE (2-4 hours)..."

$BIGSCAPE \
    --input-dir "${RES_BIGSCAPE}/input" \
    --output-dir "${RES_BIGSCAPE}/raw_output" \
    --pfam-path "$PFAM_HMM" \
    --mibig-version 3.1 \
    --cutoffs 0.3 0.5 0.7 \
    --include-singletons \
    --cpus "$N_CPUS" \
    --verbose

echo "[09.3] BiG-SCAPE complete."
find "${RES_BIGSCAPE}/raw_output" -type f | wc -l | xargs echo "  Output files:"
find "${RES_BIGSCAPE}/raw_output" -name "*.tsv" | wc -l | xargs echo "  TSV files:"

echo "=== [09] Done: $(date) ==="
