#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/13_make_figures.R
# ============================================================
# Purpose: Generate publication-grade figures for the manuscript.
#
# Outputs (into figures/):
#   Fig2_genome_landscape.{pdf,png}
#   Fig3_bgc_classes.{pdf,png}
#   Fig4_novelty_landscape.{pdf,png}
#   Fig5_phylogenetic_structure.{pdf,png}
#   Fig6_gcf_analysis.{pdf,png}
#   Fig7_ecological_phylostrat.{pdf,png}
#   Fig8_prioritized_candidates.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(ape); library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(cowplot); library(viridis); library(RColorBrewer)
  library(scales)
})

PALETTE <- list(
  primary   = "#2E86AB",
  secondary = "#A23B72",
  accent    = "#F18F01",
  rare      = "#C73E1D",
  core      = "#3B8EA5",
  neutral   = "#333333"
)

theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
  theme(
    plot.title = element_text(face = "bold", size = base_size + 1),
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = base_size - 1),
    legend.title = element_text(size = base_size - 1, face = "bold"),
    legend.text = element_text(size = base_size - 2),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(6, 6, 6, 6)
  )
}

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(root, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
regions <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_regions.tsv"), stringsAsFactors = FALSE)
classes <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv"), stringsAsFactors = FALSE)
novelty <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv"), stringsAsFactors = FALSE)

# ============================================================
# Figure 2 — Genome landscape
# ============================================================
cat("Fig 2: Genome landscape...\n")

bgc_counts <- regions %>%
  group_by(genome_accession) %>%
  summarise(n_bgcs = n(), .groups = "drop") %>%
  arrange(desc(n_bgcs))
bgc_counts$rank <- seq_len(nrow(bgc_counts))

pA <- ggplot(bgc_counts, aes(rank, n_bgcs)) +
  geom_col(fill = PALETTE$primary, width = 1) +
  geom_hline(yintercept = mean(bgc_counts$n_bgcs), linetype = "dashed", color = PALETTE$accent) +
  labs(x = "Genome (sorted)", y = "BGC regions per genome",
       title = "A  BGC abundance across genomes") + theme_pub()

pB <- ggplot(bgc_counts, aes(n_bgcs)) +
  geom_histogram(binwidth = 1, fill = PALETTE$secondary, color = "white") +
  labs(x = "BGC regions per genome", y = "Genomes",
       title = "B  Distribution") + theme_pub()

fig2 <- pA / pB + plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "Figure 2. BGC abundance across Devosia genomes",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(fig_dir, "Fig2_genome_landscape.pdf"), fig2, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_genome_landscape.png"), fig2, width = 9, height = 6, dpi = 300)

# ============================================================
# Figure 3 — BGC class composition
# ============================================================
cat("Fig 3: BGC class composition...\n")

class_counts <- classes %>% count(product_class, sort = TRUE)

p3A <- ggplot(class_counts, aes(reorder(product_class, n), n)) +
  geom_col(fill = PALETTE$primary) +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(x = "", y = "Number of class assignments",
       title = "A  Biosynthetic class distribution") +
  theme_pub() +
  expand_limits(y = max(class_counts$n) * 1.15)

# Genome x class heatmap
genome_class <- classes %>%
  count(genome_accession, product_class)

p3B <- ggplot(genome_class, aes(product_class, genome_accession, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Count") +
  labs(x = "", y = "", title = "B  Genome × class heatmap") +
  theme_pub(base_size = 7) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)
  )

fig3 <- p3A / p3B + plot_layout(heights = c(1, 2)) +
  plot_annotation(
    title = "Figure 3. Biosynthetic class composition",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(fig_dir, "Fig3_bgc_classes.pdf"), fig3, width = 10, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "Fig3_bgc_classes.png"), fig3, width = 10, height = 10, dpi = 300)

# ============================================================
# Figure 4 — Novelty landscape
# ============================================================
cat("Fig 4: Novelty landscape...\n")

novelty_summary <- novelty %>% count(novelty_category) %>%
  mutate(pct = 100 * n / sum(n))

p4A <- ggplot(novelty_summary, aes(reorder(novelty_category, -n), n, fill = novelty_category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(
    putatively_novel = PALETTE$accent,
    related = PALETTE$secondary,
    known = PALETTE$primary
  )) +
  labs(x = "", y = "BGCs", title = "A  Novelty classification") +
  theme_pub() +
  expand_limits(y = max(novelty_summary$n) * 1.15)

with_matches <- novelty %>% filter(n_kcb_hits > 0)
p4B <- ggplot(with_matches, aes(max_kcb_similarity)) +
  geom_histogram(bins = 30, fill = PALETTE$primary, color = "white") +
  geom_vline(xintercept = 70, linetype = "dashed", color = PALETTE$accent, linewidth = 0.7) +
  labs(x = "KnownClusterBlast similarity (%)",
       y = "BGCs",
       title = "B  Similarity distribution") +
  theme_pub()

fig4 <- p4A | p4B +
  plot_annotation(
    title = "Figure 4. BGC novelty landscape",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(fig_dir, "Fig4_novelty_landscape.pdf"), fig4, width = 11, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig4_novelty_landscape.png"), fig4, width = 11, height = 5, dpi = 300)

# ============================================================
# Figure 5 — Phylogenetic structure
# ============================================================
cat("Fig 5: Phylogenetic structure...\n")

lambda_file <- file.path(root, "results/11_stats/pagels_lambda.tsv")
logreg_file <- file.path(root, "results/11_stats/logistic_regression.tsv")

if (file.exists(lambda_file)) {
  lambda_df <- read.delim(lambda_file, stringsAsFactors = FALSE)

  p5A <- ggplot(lambda_df, aes(trait, lambda)) +
    geom_col(fill = PALETTE$primary, width = 0.6) +
    geom_text(aes(label = sprintf("%.3f\nP = %s", lambda, p_value)),
              vjust = -0.3, size = 3) +
    ylim(0, 1.15) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    labs(x = "", y = "Pagel's lambda",
         title = "A  Phylogenetic signal") +
    theme_pub()

  if (file.exists(logreg_file)) {
    logreg <- read.delim(logreg_file, stringsAsFactors = FALSE)
    sig <- logreg %>% filter(fdr < 0.05) %>% arrange(beta)

    if (nrow(sig) > 0) {
      p5B <- ggplot(sig, aes(beta, reorder(class, beta))) +
        geom_point(size = 3, color = PALETTE$secondary) +
        geom_errorbarh(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se),
                       height = 0.2, color = PALETTE$secondary) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        labs(x = "Logistic regression coefficient (β)",
             y = "",
             title = "B  Significant class associations (FDR<0.05)") +
        theme_pub()

      fig5 <- p5A | p5B +
        plot_annotation(
          title = "Figure 5. Phylogenetic structure of BGC repertoire",
          theme = theme(plot.title = element_text(face = "bold", size = 13))
        )
      ggsave(file.path(fig_dir, "Fig5_phylogenetic_structure.pdf"), fig5, width = 12, height = 5, dpi = 300)
      ggsave(file.path(fig_dir, "Fig5_phylogenetic_structure.png"), fig5, width = 12, height = 5, dpi = 300)
    }
  }
}

# ============================================================
# Figure 6 — GCF analysis
# ============================================================
cat("Fig 6: GCF analysis...\n")

gcf_file <- file.path(root, "results/09_bigscape/parsed/gcf_prevalence_c0.7.tsv")
if (file.exists(gcf_file)) {
  gcf <- read.delim(gcf_file, stringsAsFactors = FALSE)

  p6A <- ggplot(gcf, aes(prevalence)) +
    geom_histogram(bins = 50, fill = PALETTE$primary, color = "white") +
    geom_vline(xintercept = c(10, 90), linetype = "dashed",
               color = c(PALETTE$rare, PALETTE$core), linewidth = 0.6) +
    labs(x = "GCF prevalence (%)", y = "Number of GCFs",
         title = "A  Prevalence distribution (c0.7)") +
    theme_pub()

  p6B <- ggplot(gcf, aes(reorder(gcf_id, prevalence), prevalence, fill = category)) +
    geom_col() +
    scale_fill_manual(values = c(core = PALETTE$core,
                                  intermediate = PALETTE$secondary,
                                  rare = PALETTE$rare),
                      name = "Category") +
    labs(x = "GCFs (sorted by prevalence)", y = "Prevalence (%)",
         title = "B  GCFs ranked by prevalence") +
    theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  fig6 <- p6A | p6B +
    plot_annotation(
      title = "Figure 6. GCF diversity and prevalence",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  ggsave(file.path(fig_dir, "Fig6_gcf_analysis.pdf"), fig6, width = 12, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "Fig6_gcf_analysis.png"), fig6, width = 12, height = 5, dpi = 300)
}

# ============================================================
# Figure 7 — Ecological + phylostratigraphy
# ============================================================
cat("Fig 7: Ecological + phylostratigraphy...\n")

eco_file <- file.path(root, "results/11_stats/ecological_association_classes.tsv")
age_file <- file.path(root, "results/11_stats/gcf_ages.tsv")

panels <- list()

if (file.exists(eco_file)) {
  eco <- read.delim(eco_file, stringsAsFactors = FALSE)
  sig_eco <- eco %>% filter(fdr < 0.05)

  if (nrow(sig_eco) > 0) {
    p7A <- ggplot(sig_eco, aes(niche, class, fill = -log10(fdr))) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "plasma", direction = -1,
                            name = "-log10(FDR)") +
      labs(x = "", y = "", title = "A  Niche × class associations") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    panels[[length(panels) + 1]] <- p7A
  }
}

if (file.exists(age_file)) {
  age <- read.delim(age_file, stringsAsFactors = FALSE)
  age_counts <- age %>% count(age_category)

  p7B <- ggplot(age_counts, aes(reorder(age_category, -n), n, fill = age_category)) +
    geom_col(show.legend = FALSE) +
    scale_fill_viridis_d(option = "viridis") +
    labs(x = "Age category", y = "Number of GCFs",
         title = "B  GCF evolutionary ages") +
    theme_pub()
  panels[[length(panels) + 1]] <- p7B
}

if (length(panels) > 0) {
  fig7 <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title = "Figure 7. Ecological and evolutionary structure",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  ggsave(file.path(fig_dir, "Fig7_ecological_phylostrat.pdf"), fig7, width = 12, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "Fig7_ecological_phylostrat.png"), fig7, width = 12, height = 5, dpi = 300)
}

# ============================================================
# Figure 8 — Prioritized candidates
# ============================================================
cat("Fig 8: Prioritized candidates...\n")

cand_file <- file.path(root, "results/12_prioritization/prioritized_candidates.tsv")
if (file.exists(cand_file)) {
  cand <- read.delim(cand_file, stringsAsFactors = FALSE)

  p8A <- ggplot(cand, aes(reorder(bgc_id, score), score, fill = novelty_category)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c(
      putatively_novel = PALETTE$accent,
      related = PALETTE$secondary,
      known = PALETTE$primary
    ), name = "Novelty") +
    labs(x = "", y = "Prioritization score",
         title = "A  Top 15 candidate BGCs") +
    theme_pub()

  fig8 <- p8A +
    plot_annotation(
      title = "Figure 8. Prioritized BGCs for experimental follow-up",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  ggsave(file.path(fig_dir, "Fig8_prioritized_candidates.pdf"), fig8, width = 10, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, "Fig8_prioritized_candidates.png"), fig8, width = 10, height = 6, dpi = 300)
}

cat("\n=== [13] All figures generated ===\n")
cat(sprintf("Output: %s\n", fig_dir))
print(list.files(fig_dir, pattern = "\\.(pdf|png)$"))
