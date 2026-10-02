#!/usr/bin/env bash
# ============================================================
# scripts/05b_build_phylogeny.sh
# ============================================================
# Purpose: Build a genome-scale phylogeny from GTDB-Tk markers.
#
# Inputs:   results/05_gtdbtk/input_genomes/  (query_*.fna files)
# Outputs:
#   results/05b_phylogeny/identify_all/
#   results/05b_phylogeny/align_all/
#   results/05b_phylogeny/devosia.tree         (raw tree, query_ prefix)
#   results/05b_phylogeny/devosia_pruned.tree  (final tree, clean accessions)
#
# Tools:    GTDB-Tk 2.7.2 + FastTree 2.2.0
# Runtime:  30-60 min
# ============================================================
set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/05b_phylogeny_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [05b] Phylogeny: $(date) ==="
mkdir -p "$RES_PHYLO"

export GTDBTK_DATA_PATH="$GTDBTK_DATA"

# ---- Step 1: Identify markers on ALL genomes ----
if [[ ! -f "${RES_PHYLO}/identify_all/identify/gtdbtk.bac120.markers_summary.tsv" ]]; then
    echo "[05b.1] Identifying markers (~5-10 min)..."
    rm -rf "${RES_PHYLO}/identify_all"
    $GTDBTK identify \
        --genome_dir "${RES_GTDBTK}/input_genomes" \
        --out_dir "${RES_PHYLO}/identify_all" \
        --extension fna \
        --cpus "$N_CPUS"
else
    echo "[05b.1] Markers already identified, skipping"
fi

# ---- Step 2: Align markers ----
if [[ ! -f "${RES_PHYLO}/align_all/align/gtdbtk.bac120.user_msa.fasta.gz" ]]; then
    echo "[05b.2] Aligning markers (~15-30 min)..."
    rm -rf "${RES_PHYLO}/align_all"
    $GTDBTK align \
        --identify_dir "${RES_PHYLO}/identify_all/identify" \
        --out_dir "${RES_PHYLO}/align_all" \
        --cpus "$N_CPUS"
else
    echo "[05b.2] Alignment already exists, skipping"
fi

# ---- Step 3: Extract MSA ----
echo "[05b.3] Extracting MSA..."
gunzip -c "${RES_PHYLO}/align_all/align/gtdbtk.bac120.user_msa.fasta.gz" \
    > "${RES_PHYLO}/devosia_alignment.fasta"

# ---- Step 4: FastTree ----
echo "[05b.4] Building tree with FastTree..."
$FASTTREE -wag -gamma "${RES_PHYLO}/devosia_alignment.fasta" \
    > "${RES_PHYLO}/devosia.tree" 2> "${RES_PHYLO}/fasttree.log"

# ---- Step 5: Strip query_ prefix ----
echo "[05b.5] Stripping query_ prefix..."
sed 's/query_//g' "${RES_PHYLO}/devosia.tree" > "${RES_PHYLO}/devosia_pruned.tree"

# ---- Step 6: Verify ----
tip_count=$(grep -o "GCF_\|GCA_" "${RES_PHYLO}/devosia_pruned.tree" | wc -l)
echo "  Tree tips: $tip_count"

echo "=== [05b] Done: $(date) ==="
