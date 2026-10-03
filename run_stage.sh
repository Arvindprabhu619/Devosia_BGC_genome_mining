#!/usr/bin/env bash
# ============================================================
# Devosia BGC Genome Mining — Single-Stage Runner
# ============================================================
# Purpose:  Run just one pipeline stage (useful for debugging or
#           re-running a failed stage without restarting the full pipeline.
#
# Usage:    bash run_stage.sh <stage>
#
# Examples:
#   bash run_stage.sh 06      # run antiSMASH only
#   bash run_stage.sh 11      # run phylogenetic stats only
#   bash run_stage.sh 05b     # build phylogeny only
# ============================================================

set -euo pipefail
cd "$(dirname "$0")"
source scripts/00_paths.sh
export PATH="$(dirname $PYTHON_BIN):$(dirname $RSCRIPT):$PATH"

STAGE="${1:-}"
if [[ -z "$STAGE" ]]; then
    echo "Usage: bash run_stage.sh <stage>"
    echo ""
    echo "Available stages:"
    echo "  01   Download genomes from NCBI"
    echo "  02   Curate metadata"
    echo "  03   Build manifest"
    echo "  04   CheckM2 quality"
    echo "  05   GTDB-Tk taxonomy"
    echo "  05b  Build phylogeny"
    echo "  06   antiSMASH BGC prediction"
    echo "  07   Extract BGCs"
    echo "  08   Classify novelty"
    echo "  09   BiG-SCAPE"
    echo "  10   GCF prevalence"
    echo "  11   Phylogenetic stats"
    echo "  12   Prioritize candidates"
    echo "  13   Generate figures"
    echo "  14   Ecological association"
    echo "  15   Phylostratigraphy"
    exit 1
fi

case "$STAGE" in
    01)  bash     scripts/01_download_genomes.sh ;;
    02)  python3  scripts/02_curate_metadata.py ;;
    03)  bash     scripts/03_build_manifest.sh ;;
    04)  bash     scripts/04_run_checkm2.sh ;;
    05)  bash     scripts/05_run_gtdbtk.sh ;;
    05b) bash     scripts/05b_build_phylogeny.sh ;;
    06)  bash     scripts/06_run_antismash.sh ;;
    07)  python3  scripts/07_extract_bgcs.py ;;
    08)  python3  scripts/08_classify_novelty.py ;;
    09)  bash     scripts/09_run_bigscape.sh ;;
    10)  python3  scripts/10_gcf_prevalence.py ;;
    11)  "$WRAPPER_BIN/run_rscript"  scripts/11_phylogenetic_stats.R ;;
    12)  python3  scripts/12_prioritize_candidates.py ;;
    13)  "$WRAPPER_BIN/run_rscript"  scripts/13_make_figures.R ;;
    14)  python3  scripts/14_ecological_association.py ;;
    15)  "$WRAPPER_BIN/run_rscript"  scripts/15_phylostratigraphy.R ;;
    *)
        echo "Unknown stage: $STAGE"
        echo "Run 'bash run_stage.sh' for a list of stages."
        exit 1
        ;;
esac
