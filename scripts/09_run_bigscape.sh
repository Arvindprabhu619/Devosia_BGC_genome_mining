#!/usr/bin/env bash
# ============================================================
# scripts/09_run_bigscape.sh
# ============================================================
# Purpose: Cluster BGCs into Gene Cluster Families with BiG-SCAPE.
#
# Inputs:   results/06_antismash/raw_output/*/*.gbk
# Outputs:
#   results/09_bigscape/raw_output/
#   results/09_bigscape/parsed/gcf_assignments_c{0.3,0.5,0.7}.tsv
#
# Tool:     BiG-SCAPE 2.0.3
# Runtime:  2-4 hr
# ============================================================
set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/09_bigscape_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [09] BiG-SCAPE: $(date) ==="
echo "Pfam: $PFAM_HMM"

if [[ ! -f "$PFAM_HMM" ]]; then
    echo "ERROR: Pfam-A.hmm not found at $PFAM_HMM"
    exit 1
fi

mkdir -p "${RES_BIGSCAPE}/input"
mkdir -p "${RES_BIGSCAPE}/raw_output"
mkdir -p "${RES_BIGSCAPE}/parsed"

# ---- Collect GBK files ----
echo "[09.1] Collecting antiSMASH GenBank files..."
n=0
for gbk in "${RES_ANTISMASH}"/raw_output/*/*.gbk; do
    [[ -f "$gbk" ]] || continue
    # Skip region files
    [[ "$gbk" == *".region"* ]] && continue
    [[ "$gbk" == *".final"* ]] && continue
    acc=$(basename "$(dirname "$gbk")")
    ln -sf "$gbk" "${RES_BIGSCAPE}/input/${acc}.gbk"
    n=$((n + 1))
done
echo "  Linked $n GenBank files"

# ---- Run BiG-SCAPE ----
echo "[09.2] Running BiG-SCAPE (2-4 hr)..."

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

echo "=== [09] Done: $(date) ==="
