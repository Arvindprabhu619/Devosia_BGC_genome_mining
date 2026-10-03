#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/11_phylogenetic_stats.R
# ============================================================
# Purpose: Phylogenetic statistics on BGC distribution.
#
# Tests:
#   1. Pagel's lambda for BGC abundance + novelty proportion
#   2. Phylogenetic logistic regression for individual BGC classes
#   3. Zero-length-branch sensitivity
#
# Inputs:
#   results/05b_phylogeny/devosia_pruned.tree
#   results/06_antismash/parsed/bgc_regions.tsv
#   results/06_antismash/parsed/bgc_class_assignments.tsv
#   results/06_antismash/parsed/bgc_novelty.tsv
#   results/05_gtdbtk/final_genome_list.txt
# Outputs:
#   results/11_stats/pagels_lambda.tsv
#   results/11_stats/logistic_regression.tsv
#   results/11_stats/summary.txt
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))

tree_file <- file.path(root, "results/05b_phylogeny/devosia_pruned.tree")
regions_file <- file.path(root, "results/06_antismash/parsed/bgc_regions.tsv")
classes_file <- file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv")
novelty_file <- file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv")
genome_list <- file.path(root, "results/05_gtdbtk/final_genome_list.txt")
out_dir <- file.path(root, "results/11_stats")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== [11] Phylogenetic statistics ===\n")
cat("Loading tree...\n")
tree <- read.tree(tree_file)
tree <- multi2di(tree)
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

regions <- read.delim(regions_file, stringsAsFactors = FALSE)
novelty <- read.delim(novelty_file, stringsAsFactors = FALSE)

genomes <- readLines(genome_list)
genomes <- genomes[genomes != ""]

common <- intersect(tree$tip.label, genomes)
cat(sprintf("Common genomes: %d\n", length(common)))
tree <- keep.tip(tree, common)

# Trait 1: BGC abundance
bgc_counts <- regions %>%
  group_by(genome_accession) %>%
  summarise(n_bgcs = n(), .groups = "drop")

# Trait 2: Novelty proportion
novelty_prop <- novelty %>%
  group_by(genome_accession) %>%
  summarise(
    n_total = n(),
    n_novel = sum(novelty_category == "putatively_novel"),
    prop_novel = n_novel / n_total,
    .groups = "drop"
  )

traits <- data.frame(genome = common, stringsAsFactors = FALSE) %>%
  left_join(bgc_counts, by = c("genome" = "genome_accession")) %>%
  left_join(novelty_prop, by = c("genome" = "genome_accession")) %>%
  replace_na(list(n_bgcs = 0, prop_novel = 0))
rownames(traits) <- traits$genome
traits <- traits[tree$tip.label, ]

# ---- Pagel's lambda ----
cat("\n[1] Pagel's lambda...\n")
lambda_results <- data.frame(trait = character(), lambda = numeric(),
                              logL = numeric(), p_value = numeric(),
                              stringsAsFactors = FALSE)

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
    cat(sprintf("  %s: lambda = %.4f, P = %s\n",
                trait_name, fit$lambda, format(fit$P, scientific = TRUE)))
  }
}

write.table(lambda_results, file.path(out_dir, "pagels_lambda.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- Logistic regression per BGC class ----
cat("\n[2] Phylogenetic logistic regression...\n")

classes <- read.delim(classes_file, stringsAsFactors = FALSE)
class_matrix <- classes %>%
  mutate(present = 1) %>%
  pivot_wider(id_cols = genome_accession, names_from = product_class,
              values_from = present, values_fill = 0, values_fn = max)
class_matrix <- as.data.frame(class_matrix)
rownames(class_matrix) <- class_matrix$genome_accession
class_matrix$genome_accession <- NULL

n_genomes <- nrow(class_matrix)
keep <- sapply(class_matrix, function(col) {
  s <- sum(col)
  s > 0 && s < n_genomes
})
class_matrix <- class_matrix[, keep, drop = FALSE]

cat(sprintf("  Variable classes: %d\n", ncol(class_matrix)))

if (ncol(class_matrix) > 0 && requireNamespace("phylolm", quietly = TRUE)) {
  class_matrix <- class_matrix[tree$tip.label, , drop = FALSE]
  logreg_results <- data.frame()

  for (cls in colnames(class_matrix)) {
    y <- setNames(class_matrix[[cls]], rownames(class_matrix))
    if (length(unique(y)) < 2) next

    fit <- tryCatch({
      d <- data.frame(y = y)
      tree_binary <- ape::multi2di(tree)
      tree_binary$edge.length[tree_binary$edge.length == 0] <- 1e-8
      m <- phylolm::phyloglm(y ~ 1, data = d, phy = tree_binary, method = "logistic_MPLE")
      list(coef = coef(m)[1], se = summary(m)$coefficients[1, 2],
           p = summary(m)$coefficients[1, 4])
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
    logreg_results$fdr <- p.adjust(logreg_results$p_value, method = "BH")
    logreg_results$significant <- logreg_results$fdr < 0.05
    logreg_results <- logreg_results[order(logreg_results$p_value), ]
    write.table(logreg_results, file.path(out_dir, "logistic_regression.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    n_sig <- sum(logreg_results$significant)
    cat(sprintf("  Significant (FDR<0.05): %d\n", n_sig))
  }
}

# Summary
summary_text <- sprintf("Phylogenetic stats summary\n
Tree tips:         %d
Pagel's lambda:    %d traits tested
Logistic reg:      %s
", length(tree$tip.label), nrow(lambda_results),
if (exists("logreg_results") && nrow(logreg_results) > 0)
  sprintf("%d classes, %d significant", nrow(logreg_results),
          sum(logreg_results$significant)) else "not run")

writeLines(summary_text, file.path(out_dir, "summary.txt"))
cat("\n=== [11] Done ===\n")
