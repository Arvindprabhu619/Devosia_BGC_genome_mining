#!/usr/bin/env bash
# ============================================================
# Devosia BGC Genome Mining — Master Pipeline
# ============================================================
# Purpose:  Run the complete genome mining pipeline end-to-end.
# Usage:    bash run_pipeline.sh [--from STAGE] [--to STAGE]
#
# Stages:
#   01   Download Devosia genomes from NCBI
#   02   Curate metadata (remove MAGs, duplicates, fragmented)
#   03   Build genome manifest
#   04   CheckM2 quality assessment
#   05   GTDB-Tk taxonomy validation
#   05b  Build phylogeny (align + FastTree)
#   06   antiSMASH BGC prediction
#   07   Extract BGCs from antiSMASH output
#   08   Classify novelty vs MIBiG
#   09   BiG-SCAPE gene cluster family analysis
#   10   GCF prevalence computation
#   11   Phylogenetic statistics (Pagel's lambda, logistic regression)
#   12   Prioritize candidate BGCs
#   13   Generate publication figures
#   14   Ecological association analysis
#   15   Phylostratigraphy (BGC evolutionary ages)
#
# Runtime:  ~30 hours on 16-core / 128 GB node
# ============================================================

set -euo pipefail
cd "$(dirname "$0")"
source scripts/00_paths.sh
export PATH="$(dirname $PYTHON_BIN):$(dirname $RSCRIPT):$PATH"

# ---------- Parse arguments ----------
FROM="${1:-01}"
TO="${2:-15}"

# ---------- Setup log ----------
mkdir -p logs
LOG="logs/run_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "Devosia BGC Genome Mining Pipeline"
echo "Started:   $(date)"
echo "Range:     stages ${FROM} → ${TO}"
echo "Log file:  ${LOG}"
echo "============================================================"
echo ""

# ---------- Helper: run one stage ----------
run_stage() {
    local stage="$1"
    local script="$2"
    local description="$3"

    echo ""
    echo "============================================================"
    echo "STAGE ${stage} — ${description}"
    echo "Script: ${script}"
    echo "Started: $(date)"
    echo "============================================================"

    local start_time
    start_time=$(date +%s)

    if ! "${script}"; then
        echo ""
        echo "❌ STAGE ${stage} FAILED at $(date)"
        echo "See log: ${LOG}"
        echo "To resume from this stage: bash run_pipeline.sh ${stage} ${TO}"
        exit 1
    fi

    local end_time
    end_time=$(date +%s)
    local elapsed=$((end_time - start_time))

    echo ""
    echo "✅ STAGE ${stage} complete in ${elapsed}s"
}

# ---------- Run pipeline ----------
[[ "$FROM" < "02" || "$FROM" == "01" ]] && [[ "$TO" > "00" ]] && \
    run_stage 01 scripts/01_download_genomes.sh "Download Devosia genomes from NCBI"

[[ "$FROM" < "03" || "$FROM" == "02" ]] && [[ "$TO" > "01" ]] && \
    run_stage 02 scripts/02_curate_metadata.py "Curate metadata (remove MAGs, duplicates)"

[[ "$FROM" < "04" || "$FROM" == "03" ]] && [[ "$TO" > "02" ]] && \
    run_stage 03 scripts/03_build_manifest.sh "Build genome manifest"

[[ "$FROM" < "05" || "$FROM" == "04" ]] && [[ "$TO" > "03" ]] && \
    run_stage 04 scripts/04_run_checkm2.sh "CheckM2 genome quality"

[[ "$FROM" < "05b" || "$FROM" == "05" ]] && [[ "$TO" > "04" ]] && \
    run_stage 05 scripts/05_run_gtdbtk.sh "GTDB-Tk taxonomy validation"

[[ "$FROM" < "06" || "$FROM" == "05b" ]] && [[ "$TO" > "05" ]] && \
    run_stage 05b scripts/05b_build_phylogeny.sh "Build genome phylogeny"

[[ "$FROM" < "07" || "$FROM" == "06" ]] && [[ "$TO" > "05b" ]] && \
    run_stage 06 scripts/06_run_antismash.sh "antiSMASH BGC prediction"

[[ "$FROM" < "08" || "$FROM" == "07" ]] && [[ "$TO" > "06" ]] && \
    run_stage 07 scripts/07_extract_bgcs.py "Extract BGCs"

[[ "$FROM" < "09" || "$FROM" == "08" ]] && [[ "$TO" > "07" ]] && \
    run_stage 08 scripts/08_classify_novelty.py "Classify novelty vs MIBiG"

[[ "$FROM" < "10" || "$FROM" == "09" ]] && [[ "$TO" > "08" ]] && \
    run_stage 09 scripts/09_run_bigscape.sh "BiG-SCAPE GCF clustering"

[[ "$FROM" < "11" || "$FROM" == "10" ]] && [[ "$TO" > "09" ]] && \
    run_stage 10 scripts/10_gcf_prevalence.py "GCF prevalence"

[[ "$FROM" < "12" || "$FROM" == "11" ]] && [[ "$TO" > "10" ]] && \
    run_stage 11 scripts/11_phylogenetic_stats.R "Phylogenetic statistics"

[[ "$FROM" < "13" || "$FROM" == "12" ]] && [[ "$TO" > "11" ]] && \
    run_stage 12 scripts/12_prioritize_candidates.py "Prioritize candidates"

[[ "$FROM" < "14" || "$FROM" == "13" ]] && [[ "$TO" > "12" ]] && \
    run_stage 13 scripts/13_make_figures.R "Generate figures"

[[ "$FROM" < "15" || "$FROM" == "14" ]] && [[ "$TO" > "13" ]] && \
    run_stage 14 scripts/14_ecological_association.py "Ecological association"

[[ "$TO" > "14" ]] && \
    run_stage 15 scripts/15_phylostratigraphy.R "Phylostratigraphy"

echo ""
echo "============================================================"
echo "🎉 PIPELINE COMPLETE"
echo "Finished: $(date)"
echo ""
echo "Results:  results/"
echo "Figures:  figures/"
echo "Tables:   tables/"
echo "Log:      ${LOG}"
echo "============================================================"
