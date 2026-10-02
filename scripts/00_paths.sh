#!/usr/bin/env bash
# scripts/00_paths.sh
# Central path and configuration file. Sourced by every pipeline script.

export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- Input data ---
export DATA_RAW="${PROJECT_ROOT}/data/raw_genomes"
export DATA_META="${PROJECT_ROOT}/data/metadata"
export DATA_DB="${PROJECT_ROOT}/data/databases"

# --- Results ---
export RES_CHECKM2="${PROJECT_ROOT}/results/03_checkm2"
export RES_GTDBTK="${PROJECT_ROOT}/results/04_gtdbtk"
export RES_PHYLO="${PROJECT_ROOT}/results/05_phylogeny"
export RES_ANTISMASH="${PROJECT_ROOT}/results/06_antismash"
export RES_BIGSCAPE="${PROJECT_ROOT}/results/07_bigscape"
export RES_STATS="${PROJECT_ROOT}/results/08_stats"
export RES_PRIORITY="${PROJECT_ROOT}/results/09_prioritization"

# --- Outputs ---
export DIR_FIGS="${PROJECT_ROOT}/figures"
export DIR_TABLES="${PROJECT_ROOT}/tables"
export DIR_MANUSCRIPT="${PROJECT_ROOT}/manuscript"
export DIR_LOGS="${PROJECT_ROOT}/logs"
export DIR_SCRIPTS="${PROJECT_ROOT}/scripts"
export DIR_ENV="${PROJECT_ROOT}/env"
export DIR_DOCS="${PROJECT_ROOT}/docs"

# --- Tool wrappers ---
export WRAPPER_BIN="${PROJECT_ROOT}/env/tool_wrappers"
export ANTI_SMASH="${WRAPPER_BIN}/run_antismash"
export CHECKM2="${WRAPPER_BIN}/run_checkm2"
export GTDBTK="${WRAPPER_BIN}/run_gtdbtk"
export BIGSCAPE="${WRAPPER_BIN}/run_bigscape"
export FASTTREE="${WRAPPER_BIN}/run_fasttree"
export DATASETS="${WRAPPER_BIN}/run_datasets"

# --- Databases ---
export GTDBTK_DATA="$HOME/miniconda3/envs/gtdbtk/share/gtdbtk-2.7.2"
export CHECKM2_DB="${DATA_DB}/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
export PFAM_HMM="${DATA_DB}/Pfam-A.hmm"
export MIBIG_JSON_DIR="$HOME/miniconda3/envs/antismash/lib/python3.10/site-packages/antismash/databases/clustercompare/mibig"

# --- Compute ---
export N_CPUS=16
export MAX_MEM_GB=120

# --- Create dirs if missing ---
mkdir -p "$DATA_RAW" "$DATA_META" "$DATA_DB" \
         "$RES_CHECKM2" "$RES_GTDBTK" "$RES_PHYLO" \
         "$RES_ANTISMASH" "$RES_BIGSCAPE" "$RES_STATS" "$RES_PRIORITY" \
         "$DIR_FIGS" "$DIR_TABLES" "$DIR_MANUSCRIPT" "$DIR_LOGS"

echo "[00_paths] PROJECT_ROOT = $PROJECT_ROOT"
