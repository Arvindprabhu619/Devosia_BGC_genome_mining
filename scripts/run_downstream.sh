#!/usr/bin/env bash
# ============================================================
# scripts/run_downstream.sh
# Master script: runs all downstream stages after antiSMASH.
# ============================================================
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/00_paths.sh

echo "============================================================"
echo "Devosia BGC - Downstream Pipeline"
echo "Started: $(date)"
echo "============================================================"

# Stage 07: Extract BGCs
echo ""
echo "[STAGE 07] Extracting BGCs from antiSMASH JSON..."
~/miniconda3/envs/checkm2/bin/python scripts/07_extract_bgcs.py 2>&1 | tail -20

# Stage 08: Classify novelty
echo ""
echo "[STAGE 08] Classifying novelty vs MIBiG..."
~/miniconda3/envs/checkm2/bin/python scripts/08_classify_novelty.py 2>&1 | tail -15

# Stage 09: BiG-SCAPE
echo ""
echo "[STAGE 09] Running BiG-SCAPE (this takes 2-4 hours)..."
bash scripts/09_run_bigscape.sh 2>&1 | tail -20

# Stage 10: GCF prevalence
echo ""
echo "[STAGE 10] Computing GCF prevalence..."
~/miniconda3/envs/checkm2/bin/python scripts/10_gcf_prevalence.py 2>&1 | tail -20

# Stage 11: Phylogenetic stats
echo ""
echo "[STAGE 11] Running phylogenetic statistics..."
$WRAPPER_BIN/run_rscript scripts/11_phylogenetic_stats.R 2>&1 | tail -20

# Stage 12: Prioritize candidates
echo ""
echo "[STAGE 12] Prioritizing candidates..."
~/miniconda3/envs/checkm2/bin/python scripts/12_prioritize_candidates.py 2>&1 | tail -10

# Stage 14: Ecological association
echo ""
echo "[STAGE 14] Ecological association..."
~/miniconda3/envs/checkm2/bin/python scripts/14_ecological_association.py 2>&1 | tail -15

# Stage 15: Phylostratigraphy
echo ""
echo "[STAGE 15] Phylostratigraphy..."
$WRAPPER_BIN/run_rscript scripts/15_phylostratigraphy.R 2>&1 | tail -15

# Stage 13: Figures
echo ""
echo "[STAGE 13] Generating publication figures..."
$WRAPPER_BIN/run_rscript scripts/13_make_figures.R 2>&1 | tail -20

echo ""
echo "============================================================"
echo "Downstream pipeline complete: $(date)"
echo "============================================================"
