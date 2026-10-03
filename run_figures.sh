#!/usr/bin/env bash
# ============================================================
# run_figures.sh — Regenerate all publication figures
# ============================================================
set -euo pipefail
cd "$(dirname "$0")"
# Detect working R
if ! command -v Rscript &> /dev/null; then
  echo "ERROR: Rscript not found. Load R or activate conda env."
  exit 1
fi
echo "Using R: $(which Rscript)"
$(which Rscript) --version | head -1
FIG_DIR="scripts/figures"
mkdir -p logs
LOG="logs/run_figures_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=============================================="
echo "Devosia BGC — Figure generation"
echo "Started: $(date)"
echo "Log:     ${LOG}"
echo "=============================================="

run_fig() {
  echo ""
  echo "=== $1 ==="
  if Rscript "$FIG_DIR/$2"; then
    echo "  ✅ done"
  else
    echo "  ⚠️  FAILED (see log)"
  fi
}

echo "[1/14] Figure 1  — Workflow (manual, BioRender)"
echo "  ⚠️  Fig 1 is created manually in BioRender — skip"
run_fig "Fig 2  — Phylogenetic landscape"    "02_phylogeny_itol.R"
run_fig "Fig 3  — Biosynthetic diversity"    "03_biosynthetic_diversity.R"
run_fig "Fig 4  — GCF diversity"             "04_gcf_diversity.R"
run_fig "Fig 5  — Phylogenetic structure"    "05_phylo_structure.R"
run_fig "Fig 6  — Ecological + evolutionary" "06_ecological_evolutionary.R"
run_fig "Fig 7  — Prioritized candidates"    "07_prioritized_candidates.R"
run_fig "Fig 8  — BGC neighborhoods"         "08_bgc_neighborhoods.R"
run_fig "S1     — Multitrack tree"           "S1_multitrack_tree.R"
run_fig "S2     — CheckM2 QC"                "S2_checkm2_qc.R"
run_fig "S3     — Candidate details"         "S3_candidate_details.R"
run_fig "S4     — Novelty by genus"          "S4_novelty_by_genus.R"
run_fig "S5     — GCF rarefaction"           "S5_rarefaction.R"
run_fig "S6     — BGC length distribution"   "S6_bgc_lengths.R"

echo ""
echo "=============================================="
echo "  ✅ All figures regenerated: figures/"
echo "  Log: ${LOG}"
echo "=============================================="
