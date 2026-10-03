#!/usr/bin/env Rscript
# =============================================================================
# 02_phylogeny_itol.R
# =============================================================================
# Fig 2 — Circular phylogeny with 20 BGC class rings.
#
# This figure is generated EXTERNALLY using iTOL (Interactive Tree Of Life).
# It cannot be reliably reproduced in R because ggtree 3.x requires
# ggplot2 >= 4.0.0, which has a known S7 + `+.gg` operator conflict
# that breaks other figure scripts in this repository.
#
# Instead, use the pre-generated iTOL inputs:
#
#   Tree:          results/itOL_exports/itol_tree.nwk
#   Datasets:      results/itOL_exports/itol_dataset_*.txt
#                  results/itOL_exports/evolview_dataset_*.txt
#
# Steps to reproduce:
#   1. Go to https://itol.embl.de/
#   2. Upload results/itOL_exports/itol_tree.nwk
#   3. Drag-and-drop the dataset files onto the tree
#   4. Set tree type → Circular
#   5. Enable legends (top-left for Genus, left for BGC classes)
#   6. Export as PDF → save as figures/Fig2_phylogeny_itol.pdf
#
# If you must regenerate in R, you need an environment with:
#   ggtree >= 3.10, ggtreeExtra, ggnewscale, ggplot2 >= 4.0.0
#   See docs/optional_fig2_ggtree_env.md
# =============================================================================

cat("Fig 2 — circular phylogeny is generated in iTOL.\n")
cat("Inputs:  results/itOL_exports/\n")
cat("Output:  figures/Fig2_phylogeny_itol.pdf\n")
cat("See:     docs/fig2_workflow.md\n")
