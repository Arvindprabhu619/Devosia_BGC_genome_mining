#!/usr/bin/env bash
# ============================================================
# scripts/09_run_bigscape.sh
# ============================================================
# Purpose: Cluster BGCs into GCFs with BiG-SCAPE 2.0.3.
#
# Input:  antiSMASH 7 per-region GBK files (one per BGC).
# Output: BiG-SCAPE cluster results (DB + TSV per class/cutoff).
#
# Runtime: ~30 min on 16-core node
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

rm -rf "${RES_BIGSCAPE}/input"
rm -rf "${RES_BIGSCAPE}/raw_output"
mkdir -p "${RES_BIGSCAPE}/input"
mkdir -p "${RES_BIGSCAPE}/raw_output"

# ---- Collect region-level GBK files ----
# Use full basename to preserve uniqueness across contigs within a genome.
echo "[09.1] Collecting antiSMASH region GBK files..."
n=0
for gbk in "${RES_ANTISMASH}"/raw_output/*/*.region*.gbk; do
    [[ -f "$gbk" ]] || continue
    acc=$(basename "$(dirname "$gbk")")
    base=$(basename "$gbk" .gbk)   # e.g., JAYRZA010000049.1.region001
    ln -sf "$gbk" "${RES_BIGSCAPE}/input/${acc}__${base}.gbk"
    n=$((n + 1))
done
echo "  Linked $n region GBK files"

if [[ $n -lt 600 ]]; then
    echo "ERROR: Expected ~663 region files, got only $n"
    exit 1
fi

# ---- Run BiG-SCAPE 2.0 ----
echo "[09.2] Running BiG-SCAPE cluster (30 min)..."

$BIGSCAPE cluster \
    --input-dir "${RES_BIGSCAPE}/input" \
    --output-dir "${RES_BIGSCAPE}/raw_output" \
    --pfam-path "$PFAM_HMM" \
    -m 3.1 \
    --gcf-cutoffs 0.3,0.5,0.7 \
    --include-singletons \
    --include-gbk "*" \
    --classify category \
    -c "$N_CPUS" \
    -v

echo "[09.3] BiG-SCAPE complete."
find "${RES_BIGSCAPE}/raw_output" -type f | wc -l | xargs echo "  Output files:"

echo "=== [09] Done: $(date) ==="
