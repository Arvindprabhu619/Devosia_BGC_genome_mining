#!/usr/bin/env Rscript
# ============================================================
# scripts/15_phylostratigraphy.R
#
# Purpose: Assign each GCF to an evolutionary age using
#          ancestral state reconstruction on the phylogeny.
#
# Inputs:
#   results/05_phylogeny/devosia_pruned.tree
#   results/07_bigscape/parsed/gcf_prevalence_c0.7.tsv
# Outputs:
#   results/08_stats/gcf_ages.tsv
#   results/08_stats/phylostratigraphy_summary.txt
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))

tree_file <- file.path(root, "results/05_phylogeny/devosia_pruned.tree")
gcf_file <- file.path(root, "results/07_bigscape/parsed/gcf_prevalence_c0.7.tsv")
out_dir <- file.path(root, "results/08_stats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== [15] Phylostratigraphy ===\n")

tree <- read.tree(tree_file)
tree <- multi2di(tree)
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

# Load GCF prevalence
gcf <- read.delim(gcf_file, stringsAsFactors = FALSE)

# For each GCF, use presence/absence across tips
# Run ancestral state reconstruction (equal rates, parsimony or ML)
age_results <- data.frame(
  gcf_id = character(),
  n_genomes = integer(),
  mrca_node = integer(),
  mrca_depth = numeric(),
  age_category = character(),
  stringsAsFactors = FALSE
)

# Depth of each node from root (in branch length units)
node_depths <- node.depth.edgelength(tree)
names(node_depths) <- seq_len(length(node_depths))

root_depth <- 0  # root is always 0

n_gcfs <- nrow(gcf)
cat(sprintf("Analyzing %d GCFs...\n", n_gcfs))

for (i in seq_len(n_gcfs)) {
  row <- gcf[i, ]
  gcf_id <- row$gcf_id
  genomes_str <- row$genomes
  present_genomes <- strsplit(genomes_str, ",")[[1]]
  present_genomes <- intersect(present_genomes, tree$tip.label)

  if (length(present_genomes) < 2) {
    age_results <- rbind(age_results, data.frame(
      gcf_id = gcf_id,
      n_genomes = length(present_genomes),
      mrca_node = NA,
      mrca_depth = NA,
      age_category = "singleton",
      stringsAsFactors = FALSE
    ))
    next
  }

  # MRCA node
  mrca <- tryCatch(getMRCA(tree, present_genomes), error = function(e) NA)
  if (is.na(mrca)) next

  depth <- node_depths[mrca]

  # Age category based on relative depth from root
  max_depth <- max(node_depths)
  relative_depth <- depth / max_depth

  if (relative_depth < 0.25) age_cat <- "ancient"
  else if (relative_depth < 0.50) age_cat <- "old"
  else if (relative_depth < 0.75) age_cat <- "intermediate"
  else age_cat <- "recent"

  age_results <- rbind(age_results, data.frame(
    gcf_id = gcf_id,
    n_genomes = length(present_genomes),
    mrca_node = mrca,
    mrca_depth = round(depth, 4),
    age_category = age_cat,
    stringsAsFactors = FALSE
  ))

  if (i %% 10 == 0) cat(sprintf("  Processed %d/%d\n", i, n_gcfs))
}

# Write output
write.table(age_results, file.path(out_dir, "gcf_ages.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Summary
age_counts <- table(age_results$age_category)
summary_lines <- c(
  "=== GCF Phylostratigraphy Summary ===",
  sprintf("Total GCFs: %d", nrow(age_results)),
  "",
  "Age category distribution:"
)
for (cat_name in names(age_counts)) {
  summary_lines <- c(summary_lines,
                     sprintf("  %-15s %d", cat_name, age_counts[cat_name]))
}

cat(paste(summary_lines, collapse = "\n"))
writeLines(summary_lines, file.path(out_dir, "phylostratigraphy_summary.txt"))

cat("\n=== [15] Done ===\n")
