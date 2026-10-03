#!/usr/bin/env Rscript
# =============================================================================
# S3_candidate_details.R — Supplementary: Candidate details table
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
root     <- base_dir
fig_dir  <- file.path(base_dir, "figures")

cand <- read.delim(file.path(root, "results/12_prioritization/prioritized_candidates.tsv"),
                   stringsAsFactors = FALSE)

# Simple table-as-figure visualization
cand$label <- paste0(cand$rank, ". ", cand$genome_accession, "  ", cand$product_classes)

p <- ggplot(cand, aes(x = score, y = reorder(label, score))) +
  geom_col(fill = "#2166AC", width = 0.7) +
  geom_text(aes(label = sprintf("%s (score %d)", novelty_category, score)),
            hjust = -0.05, size = 3) +
  scale_x_continuous(limits = c(0, 8)) +
  labs(x = "Score", y = NULL,
       title = "Prioritized candidate BGCs") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"))

ggsave(file.path(fig_dir, "S3_candidate_details.pdf"), p,
       width = 10, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "S3_candidate_details.png"), p,
       width = 10, height = 5, dpi = 300)
cat("S3 saved\n")
