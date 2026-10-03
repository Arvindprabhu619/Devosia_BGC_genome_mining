#!/usr/bin/env Rscript
# =============================================================================
# 25_make_rarefaction.R
# Supplementary: GCF rarefaction curve
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
in_file  <- file.path(base_dir, "results/09_bigscape/parsed/gcf_prevalence_c0.7.tsv")
out_dir  <- file.path(base_dir, "figures/Supplementary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gcf <- read.delim(in_file, stringsAsFactors = FALSE)
cat("GCFs:", nrow(gcf), "\n")

# Parse comma-separated genome lists into long format
genomes_long <- gcf %>%
  select(gcf_id, genomes) %>%
  separate_rows(genomes, sep = ",") %>%
  mutate(genomes = trimws(genomes)) %>%
  filter(genomes != "") %>%
  distinct()

cat("Total genome-GCF pairs:", nrow(genomes_long), "\n")
cat("Unique genomes:", length(unique(genomes_long$genomes)), "\n")

# Build presence matrix: rows = GCFs, cols = genomes
presence <- genomes_long %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = genomes, values_from = present, values_fill = 0) %>%
  as.data.frame()
rownames(presence) <- presence$gcf_id
presence <- presence[, -1, drop = FALSE]

cat("Presence matrix:", nrow(presence), "x", ncol(presence), "\n")

# Rarefaction
set.seed(42)
n_genomes <- ncol(presence)
n_iter    <- 100

rarefy <- function(n) {
  cols <- sample(n_genomes, n)
  sum(rowSums(presence[, cols, drop = FALSE]) > 0)
}

rar_df <- expand.grid(n_genomes = 1:n_genomes, iter = 1:n_iter) %>%
  rowwise() %>%
  mutate(n_gcfs = rarefy(n_genomes)) %>%
  ungroup() %>%
  group_by(n_genomes) %>%
  summarise(
    mean_gcfs = mean(n_gcfs),
    lo        = quantile(n_gcfs, 0.025),
    hi        = quantile(n_gcfs, 0.975),
    .groups   = "drop"
  )

p <- ggplot(rar_df, aes(x = n_genomes, y = mean_gcfs)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#67A9CF", alpha = 0.4) +
  geom_line(color = "#2166AC", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, n_genomes, by = 20)) +
  labs(
    title = "GCF rarefaction curve (c0.7)",
    x = "Number of genomes sampled",
    y = "Cumulative unique GCFs",
    caption = "Mean + 95% CI from 100 random samplings"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "S_rarefaction_GCF.pdf"), p, width = 8, height = 5)
ggsave(file.path(out_dir, "S_rarefaction_GCF.png"), p, width = 8, height = 5, dpi = 300)
cat("Saved:", out_dir, "\n")
