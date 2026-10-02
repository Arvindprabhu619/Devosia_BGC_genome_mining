#!/usr/bin/env bash
# ============================================================
# scripts/00_paths.sh
# ============================================================
# Purpose:
#   Central path and configuration file. Sourced by every pipeline
#   script to ensure consistent directory structure and tool paths.
#
# Usage:
#   source "$(dirname "$0")/00_paths.sh"
#
# Notes:
#   - All paths are relative to PROJECT_ROOT, auto-detected
#   - Creates all necessary output directories if missing
#   - Exports tool wrapper paths
#   - Exports database paths
#   - Exports GTDBTK_DATA_PATH (required by GTDB-Tk)
# ============================================================

# ---------- Project root ----------
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---------- Input data ----------
export DATA_RAW="${PROJECT_ROOT}/data/raw_genomes"
export DATA_META="${PROJECT_ROOT}/data/metadata"
export DATA_DB="${PROJECT_ROOT}/data/databases"

# ---------- Results by stage ----------
export RES_CHECKM2="${PROJECT_ROOT}/results/04_checkm2"
export RES_GTDBTK="${PROJECT_ROOT}/results/05_gtdbtk"
export RES_PHYLO="${PROJECT_ROOT}/results/05b_phylogeny"
export RES_ANTISMASH="${PROJECT_ROOT}/results/06_antismash"
export RES_BIGSCAPE="${PROJECT_ROOT}/results/09_bigscape"
export RES_STATS="${PROJECT_ROOT}/results/11_stats"
export RES_PRIORITY="${PROJECT_ROOT}/results/12_prioritization"

# ---------- Outputs ----------
export DIR_FIGS="${PROJECT_ROOT}/figures"
export DIR_TABLES="${PROJECT_ROOT}/tables"
export DIR_MANUSCRIPT="${PROJECT_ROOT}/manuscript"
export DIR_LOGS="${PROJECT_ROOT}/logs"
export DIR_DOCS="${PROJECT_ROOT}/docs"

# ---------- Scripts ----------
export DIR_SCRIPTS="${PROJECT_ROOT}/scripts"
export DIR_UTILS="${PROJECT_ROOT}/scripts/utils"

# ---------- Tool wrappers ----------
export WRAPPER_BIN="${PROJECT_ROOT}/env/tool_wrappers"
export ANTI_SMASH="${WRAPPER_BIN}/run_antismash"
export CHECKM2="${WRAPPER_BIN}/run_checkm2"
export GTDBTK="${WRAPPER_BIN}/run_gtdbtk"
export BIGSCAPE="${WRAPPER_BIN}/run_bigscape"
export FASTTREE="${WRAPPER_BIN}/run_fasttree"
export DATASETS="${WRAPPER_BIN}/run_datasets"
export RSCRIPT="${WRAPPER_BIN}/run_rscript"
export R="${WRAPPER_BIN}/run_r"

# ---------- Python (working interpreter) ----------
# The base Python is broken; use the checkm2 env's Python
export PYTHON_BIN="${PYTHON_BIN:-$HOME/miniconda3/envs/checkm2/bin/python}"

# ---------- Databases ----------
export GTDBTK_DATA="${DATA_DB}/gtdbtk_r232/release232"
export GTDBTK_DATA_PATH="${GTDBTK_DATA}"   # required by GTDB-Tk
export CHECKM2_DB="${DATA_DB}/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
export PFAM_HMM="${DATA_DB}/Pfam-A.hmm"
export MIBIG_DIR="${DATA_DB}/mibig_3.1"

# ---------- Compute ----------
export N_CPUS="${N_CPUS:-16}"
export MAX_MEM_GB="${MAX_MEM_GB:-120}"

# ---------- Create directories if missing ----------
mkdir -p "$DATA_RAW" "$DATA_META" "$DATA_DB"
mkdir -p "$RES_CHECKM2" "$RES_GTDBTK" "$RES_PHYLO"
mkdir -p "$RES_ANTISMASH" "$RES_BIGSCAPE" "$RES_STATS" "$RES_PRIORITY"
mkdir -p "$DIR_FIGS" "$DIR_TABLES" "$DIR_MANUSCRIPT" "$DIR_LOGS"
mkdir -p "$DIR_DOCS" "$DIR_SCRIPTS" "$DIR_UTILS"

# ---------- Report ----------
if [[ "${QUIET:-0}" != "1" ]]; then
    echo "[00_paths] PROJECT_ROOT = $PROJECT_ROOT"
    echo "[00_paths] GTDBTK_DATA  = $GTDBTK_DATA"
    echo "[00_paths] CHECKM2_DB   = $CHECKM2_DB"
fi
