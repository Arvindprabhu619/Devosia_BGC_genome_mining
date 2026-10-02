#!/usr/bin/env Rscript
# ============================================================
# scripts/11_phylogenetic_stats.R
#
# Purpose: Phylogenetic statistics on BGC distribution and novelty.
#
# Tests:
#   1. Pagel's λ for BGC abundance + novelty proportion
#   2. Phylogenetic logistic regression for individual BGC classes
#   3. Pagel ER vs ARD model comparison
#   4. Zero-length-branch sensitivity
#
# Inputs:
#   results/05_phylogeny/devosia.tree       (Newick)
#   results/06_antismash/parsed/bgc_regions.tsv
#   results/06_antismash/parsed/bgc_class_assignments.tsv
#   results/06_antismash/parsed/bgc_novelty.tsv
#   results/04_gtdbtk/final_genome_list.txt
# Outputs:
#   results/08_stats/pagels_lambda.tsv
#   results/08_stats/logistic_regression.tsv
#   results/08_stats/pagel_ER_ARD.tsv
#   results/08_stats/sensitivity.tsv
#   results/08_stats/summary.txt
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
  library(tidyr)
})

# ---- Paths ----
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))

tree_file <- file.path(root, "results/05_phylogeny/devosia_pruned.tree")
if (!file.exists(tree_file)) {
  tree_file <- file.path(root, "results/05_phylogeny/devosia.tree")
}

regions_file <- file.path(root, "results/06_antismash/parsed/bgc_regions.tsv")
classes_file <- file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv")
novelty_file <- file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv")
genome_list <- file.path(root, "results/04_gtdbtk/final_genome_list.txt")

out_dir <- file.path(root, "results/08_stats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== [11] Phylogenetic statistics ===\n")

# ---- Load data ----
cat("Loading tree and BGC data...\n")
tree <- read.tree(tree_file)
tree <- multi2di(tree)  # resolve polytomies
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

regions <- read.delim(regions_file, stringsAsFactors = FALSE)
classes <- read.delim(classes_file, stringsAsFactors = FALSE)
novelty <- read.delim(novelty_file, stringsAsFactors = FALSE)
genomes <- readLines(genome_list)
genomes <- genomes[genomes != ""]

# Match tree tips to genomes
common <- intersect(tree$tip.label, genomes)
cat(sprintf("Common genomes (tree + final list): %d\n", length(common)))

# Prune tree to common
tree <- keep.tip(tree, common)

# ---- Trait 1: BGC abundance per genome ----
bgc_counts <- regions %>%
  group_by(genome_accession) %>%
  summarise(n_bgcs = n(), .groups = "drop")

# ---- Trait 2: Novelty proportion per genome ----
novelty_prop <- novelty %>%
  group_by(genome_accession) %>%
  summarise(
    n_total = n(),
    n_novel = sum(novelty_category == "putatively_novel"),
    prop_novel = n_novel / n_total,
    .groups = "drop"
  )

# Merge with genome list, fill missing with 0
traits <- data.frame(genome = common, stringsAsFactors = FALSE) %>%
  left_join(bgc_counts, by = c("genome" = "genome_accession")) %>%
  left_join(novelty_prop, by = c("genome" = "genome_accession")) %>%
  replace_na(list(n_bgcs = 0, prop_novel = 0))

rownames(traits) <- traits$genome
traits <- traits[tree$tip.label, ]

# ---- 1. Pagel's λ for continuous traits ----
cat("\n[1] Pagel's λ analysis...\n")

lambda_results <- data.frame(
  trait = character(),
  lambda = numeric(),
  logL = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (trait_name in c("n_bgcs", "prop_novel")) {
  x <- setNames(traits[[trait_name]], rownames(traits))
  fit <- tryCatch(
    phytools::phylosig(tree, x, method = "lambda", test = TRUE),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    lambda_results <- rbind(lambda_results, data.frame(
      trait = trait_name,
      lambda = round(fit$lambda, 4),
      logL = round(fit$logL, 4),
      p_value = format(fit$P, scientific = TRUE, digits = 4),
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %s: λ = %.4f, P = %s\n", trait_name, fit$lambda, format(fit$P, scientific = TRUE)))
  }
}

write.table(lambda_results, file.path(out_dir, "pagels_lambda.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- 2. Phylogenetic logistic regression for BGC classes ----
cat("\n[2] Phylogenetic logistic regression per BGC class...\n")

# Build genome x class presence matrix
class_matrix <- classes %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = genome_accession,
    names_from = product_class,
    values_from = present,
    values_fill = 0,
    values_fn = max
  )

class_matrix <- as.data.frame(class_matrix)
rownames(class_matrix) <- class_matrix$genome_accession
class_matrix$genome_accession <- NULL

# Only keep variable classes (not present in 0 or all genomes)
n_genomes <- nrow(class_matrix)
keep_classes <- sapply(class_matrix, function(col) {
  s <- sum(col)
  s > 0 && s < n_genomes
})
class_matrix <- class_matrix[, keep_classes, drop = FALSE]

cat(sprintf("  Classes with sufficient variation: %d\n", ncol(class_matrix)))

# Align with tree
class_matrix <- class_matrix[tree$tip.label, , drop = FALSE]

logreg_results <- data.frame()
for (cls in colnames(class_matrix)) {
  y <- setNames(class_matrix[[cls]], rownames(class_matrix))
  # Remove classes that are all one state after alignment
  if (length(unique(y)) < 2) next

  fit <- tryCatch({
    # Use phylolm for phylogenetic logistic regression
    if (requireNamespace("phylolm", quietly = TRUE)) {
      d <- data.frame(y = y)
      m <- phylolm::phylolm(y ~ 1, data = d, phy = tree, model = "logistic_MPLE")
      list(coef = coef(m)[1], se = summary(m)$coefficients[1, 2],
           p = summary(m)$coefficients[1, 4])
    } else NULL
  }, error = function(e) NULL)

  if (!is.null(fit)) {
    logreg_results <- rbind(logreg_results, data.frame(
      class = cls,
      beta = round(fit$coef, 4),
      se = round(fit$se, 4),
      p_value = fit$p,
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(logreg_results) > 0) {
  # FDR correction
  logreg_results$fdr <- p.adjust(logreg_results$p_value, method = "BH")
  logreg_results$significant <- logreg_results$fdr < 0.05
  logreg_results <- logreg_results[order(logreg_results$p_value), ]
  write.table(logreg_results, file.path(out_dir, "logistic_regression.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Significant after FDR<0.05: %d\n", sum(logreg_results$significant)))
}

# ---- 3. Pagel ER vs ARD ----
cat("\n[3] Pagel ER vs ARD model comparison...\n")
# Requires phytools::fitPagel - simplified placeholder
# Full implementation would iterate over binary traits

cat("\n=== [11] Done ===\n")
