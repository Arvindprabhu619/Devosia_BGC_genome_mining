#!/usr/bin/env bash
# ============================================================
# scripts/09_run_bigscape.sh
# ============================================================
# Purpose: Cluster BGCs into Gene Cluster Families with BiG-SCAPE 2.0.3
#
# Uses antiSMASH 7 per-region GBK files (region*.gbk).
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
echo "[09.1] Collecting antiSMASH region GBK files..."
n=0
for gbk in "${RES_ANTISMASH}"/raw_output/*/*.region*.gbk; do
    [[ -f "$gbk" ]] || continue
    base=$(basename "$gbk")
    acc=$(basename "$(dirname "$gbk")")
    region=$(echo "$base" | sed 's/.*\.region/region/;s/\.gbk//')
    ln -sf "$gbk" "${RES_BIGSCAPE}/input/${acc}_${region}.gbk"
    n=$((n + 1))
done
echo "  Linked $n region GBK files"

if [[ $n -eq 0 ]]; then
    echo "ERROR: No region GBK files found."
    exit 1
fi

# ---- Run BiG-SCAPE 2.0 ----
echo "[09.2] Running BiG-SCAPE cluster (2-4 hr)..."

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
