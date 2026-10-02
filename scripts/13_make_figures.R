#!/usr/bin/env Rscript
# 13_make_figures.R — Publication-grade figures for Devosia BGC paper.
suppressPackageStartupMessages({
  library(ape); library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(cowplot); library(viridis); library(RColorBrewer)
})

PALETTE <- list(
  primary = "#2E86AB", secondary = "#A23B72", accent = "#F18F01",
  rare = "#C73E1D", core = "#3B8EA5", neutral = "#333333"
)

theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
  theme(
    plot.title = element_text(face = "bold", size = base_size + 2),
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = base_size - 1),
    legend.title = element_text(size = base_size - 1, face = "bold"),
    legend.text = element_text(size = base_size - 2),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(8, 8, 8, 8)
  )
}

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(root, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")
regions <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_regions.tsv"), stringsAsFactors = FALSE)
classes <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv"), stringsAsFactors = FALSE)
novelty <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv"), stringsAsFactors = FALSE)

# Figure 2 — Genome landscape
cat("Fig 2: Genome landscape...\n")
bgc_counts <- regions %>% group_by(genome_accession) %>% summarise(n_bgcs = n(), .groups = "drop") %>% arrange(desc(n_bgcs))
bgc_counts$rank <- seq_len(nrow(bgc_counts))

pA <- ggplot(bgc_counts, aes(rank, n_bgcs)) +
  geom_col(fill = PALETTE$primary, width = 1) +
  geom_hline(yintercept = mean(bgc_counts$n_bgcs), linetype = "dashed", color = PALETTE$accent) +
  labs(x = "Genome (sorted)", y = "BGCs per genome", title = "A. BGC abundance") + theme_pub()

pB <- ggplot(bgc_counts, aes(n_bgcs)) +
  geom_histogram(binwidth = 1, fill = PALETTE$secondary, color = "white") +
  labs(x = "BGCs per genome", y = "Genomes", title = "B. Distribution") + theme_pub()

fig2 <- pA / pB + plot_layout(heights = c(2, 1)) +
  plot_annotation(title = "Figure 2. BGC abundance across 110 Devosia genomes")

ggsave(file.path(fig_dir, "Fig2_genome_landscape.pdf"), fig2, width = 10, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_genome_landscape.png"), fig2, width = 10, height = 7, dpi = 300)

# Figure 3 — Class composition
cat("Fig 3: Class composition...\n")
class_counts <- classes %>% count(product_class, sort = TRUE)
p3A <- ggplot(class_counts, aes(reorder(product_class, n), n)) +
  geom_col(fill = PALETTE$primary) +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  coord_flip() + labs(x = "", y = "Class assignments", title = "A. Biosynthetic classes") +
  theme_pub() + expand_limits(y = max(class_counts$n) * 1.1)

fig3 <- p3A + plot_annotation(title = "Figure 3. Biosynthetic class composition")
ggsave(file.path(fig_dir, "Fig3_bgc_classes.pdf"), fig3, width = 8, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "Fig3_bgc_classes.png"), fig3, width = 8, height = 7, dpi = 300)

# Figure 4 — Novelty
cat("Fig 4: Novelty landscape...\n")
novelty_summary <- novelty %>% count(novelty_category) %>% mutate(pct = 100 * n / sum(n))
p4A <- ggplot(novelty_summary, aes(reorder(novelty_category, -n), n, fill = novelty_category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(putatively_novel = PALETTE$accent, related = PALETTE$secondary, known = PALETTE$primary)) +
  labs(x = "", y = "BGCs", title = "A. Novelty classification") + theme_pub()

with_matches <- novelty %>% filter(n_kcb_hits > 0)
p4B <- ggplot(with_matches, aes(max_kcb_similarity)) +
  geom_histogram(bins = 30, fill = PALETTE$primary, color = "white") +
  geom_vline(xintercept = 70, linetype = "dashed", color = PALETTE$accent) +
  labs(x = "KCB similarity (%)", y = "BGCs", title = "B. Similarity distribution") + theme_pub()

fig4 <- p4A | p4B + plot_annotation(title = "Figure 4. BGC novelty landscape")
ggsave(file.path(fig_dir, "Fig4_novelty_landscape.pdf"), fig4, width = 11, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig4_novelty_landscape.png"), fig4, width = 11, height = 5, dpi = 300)

# Figure 6 — GCF analysis
cat("Fig 6: GCF analysis...\n")
gcf_file <- file.path(root, "results/07_bigscape/parsed/gcf_prevalence_c0.7.tsv")
if (file.exists(gcf_file)) {
  gcf <- read.delim(gcf_file, stringsAsFactors = FALSE)
  p6A <- ggplot(gcf, aes(prevalence)) +
    geom_histogram(bins = 50, fill = PALETTE$primary, color = "white") +
    geom_vline(xintercept = c(10, 90), linetype = "dashed", color = c(PALETTE$rare, PALETTE$core)) +
    labs(x = "GCF prevalence (%)", y = "GCFs", title = "A. Prevalence distribution") + theme_pub()

  p6B <- ggplot(gcf, aes(reorder(gcf_id, prevalence), prevalence, fill = category)) +
    geom_col() +
    scale_fill_manual(values = c(core = PALETTE$core, intermediate = PALETTE$secondary, rare = PALETTE$rare)) +
    labs(x = "GCFs (sorted)", y = "Prevalence (%)", title = "B. GCFs by prevalence") +
    theme_pub() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  fig6 <- p6A | p6B + plot_annotation(title = "Figure 6. GCF diversity and prevalence")
  ggsave(file.path(fig_dir, "Fig6_gcf_analysis.pdf"), fig6, width = 12, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "Fig6_gcf_analysis.png"), fig6, width = 12, height = 5, dpi = 300)
}

# Figure 7 — Prioritized candidates
cat("Fig 7: Prioritized candidates...\n")
cand_file <- file.path(root, "results/09_prioritization/prioritized_candidates.tsv")
if (file.exists(cand_file)) {
  cand <- read.delim(cand_file, stringsAsFactors = FALSE)
  p7A <- ggplot(cand, aes(reorder(bgc_id, score), score, fill = novelty_category)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c(putatively_novel = PALETTE$accent, related = PALETTE$secondary, known = PALETTE$primary)) +
    labs(x = "", y = "Score", title = "A. Top 15 candidate BGCs") + theme_pub()

  fig7 <- p7A + plot_annotation(title = "Figure 7. Prioritized candidates")
  ggsave(file.path(fig_dir, "Fig7_prioritized_candidates.pdf"), fig7, width = 9, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, "Fig7_prioritized_candidates.png"), fig7, width = 9, height = 6, dpi = 300)
}

cat("\n=== All figures generated in", fig_dir, "===\n")
print(list.files(fig_dir, pattern = "\\.(pdf|png)$"))
