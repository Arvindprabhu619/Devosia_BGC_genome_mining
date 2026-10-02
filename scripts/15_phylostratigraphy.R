#!/usr/bin/env Rscript
# ============================================================
# scripts/15_phylostratigraphy.R
# ============================================================
# Purpose: Assign GCFs to evolutionary ages via MRCA depth.
#
# Inputs:
#   results/05b_phylogeny/devosia_pruned.tree
#   results/09_bigscape/parsed/gcf_prevalence_c0.7.tsv
# Outputs:
#   results/11_stats/gcf_ages.tsv
#   results/11_stats/phylostratigraphy_summary.txt
# ============================================================

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))

tree_file <- file.path(root, "results/05b_phylogeny/devosia_pruned.tree")
gcf_file <- file.path(root, "results/09_bigscape/parsed/gcf_prevalence_c0.7.tsv")
out_dir <- file.path(root, "results/11_stats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== [15] Phylostratigraphy ===\n")

tree <- read.tree(tree_file)
tree <- multi2di(tree)
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

gcf <- read.delim(gcf_file, stringsAsFactors = FALSE)
node_depths <- node.depth.edgelength(tree)
max_depth <- max(node_depths)

age_results <- data.frame()

for (i in seq_len(nrow(gcf))) {
  row <- gcf[i, ]
  genomes_str <- row$genomes
  present <- strsplit(genomes_str, ",")[[1]]
  present <- intersect(present, tree$tip.label)

  if (length(present) < 2) {
    age_results <- rbind(age_results, data.frame(
      gcf_id = row$gcf_id, n_genomes = length(present),
      mrca_node = NA, mrca_depth = NA, age_category = "singleton",
      stringsAsFactors = FALSE
    ))
    next
  }

  mrca <- tryCatch(getMRCA(tree, present), error = function(e) NA)
  if (is.na(mrca)) next

  depth <- node_depths[mrca]
  rel <- depth / max_depth

  if (rel < 0.25) age_cat <- "ancient"
  else if (rel < 0.5) age_cat <- "old"
  else if (rel < 0.75) age_cat <- "intermediate"
  else age_cat <- "recent"

  age_results <- rbind(age_results, data.frame(
    gcf_id = row$gcf_id, n_genomes = length(present),
    mrca_node = mrca, mrca_depth = round(depth, 4),
    age_category = age_cat, stringsAsFactors = FALSE
  ))
}

write.table(age_results, file.path(out_dir, "gcf_ages.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

age_counts <- table(age_results$age_category)
summary_text <- c(
  "=== GCF Phylostratigraphy Summary ===",
  sprintf("Total GCFs: %d", nrow(age_results)),
  "",
  "Age category distribution:"
)
for (cat_name in names(age_counts)) {
  summary_text <- c(summary_text,
                    sprintf("  %-15s %d", cat_name, age_counts[cat_name]))
}
writeLines(summary_text, file.path(out_dir, "phylostratigraphy_summary.txt"))
cat(paste(summary_text, collapse = "\n"))
cat("\n=== [15] Done ===\n")
