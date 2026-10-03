#!/usr/bin/env Rscript
# =============================================================================
# 03_biosynthetic_diversity.R
# Figure 3: Biosynthetic diversity (3 panels)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
root     <- base_dir
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Palette + theme (self-contained)
PAL <- list(
  primary   = "#2166AC",
  secondary = "#B2182B",
  accent    = "#EF8A62",
  core      = "#2166AC",
  rare      = "#B2182B",
  grid      = "grey90"
)
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      strip.text       = element_text(face = "bold")
    )
}

# --- Load data ---------------------------------------------------------------
classes <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv"),
                      stringsAsFactors = FALSE)
novelty <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv"),
                      stringsAsFactors = FALSE)
genomes <- unique(classes$genome_accession)

# --- Panel A: Class composition ----------------------------------------------
class_counts <- classes %>% count(product_class, sort = TRUE) %>%
  mutate(product_class = factor(product_class, levels = rev(product_class)))

p3a <- ggplot(class_counts, aes(x = n, y = product_class)) +
  geom_col(fill = PAL$primary, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.3, size = 3) +
  labs(x = "Number of class assignments", y = NULL,
       title = "Biosynthetic class composition") +
  theme_pub() +
  xlim(0, max(class_counts$n) * 1.15)

# --- Panel B: Genome x class heatmap -----------------------------------------
genome_class <- classes %>%
  count(genome_accession, product_class) %>%
  complete(genome_accession = genomes,
           product_class = unique(classes$product_class),
           fill = list(n = 0))

gm_order <- genome_class %>%
  group_by(genome_accession) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(genome_accession)

cls_order <- levels(class_counts$product_class)
genome_class$genome_accession <- factor(genome_class$genome_accession, levels = gm_order)
genome_class$product_class <- factor(genome_class$product_class, levels = cls_order)

p3b <- ggplot(genome_class, aes(x = product_class, y = genome_accession, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Count",
                       breaks = c(1, 3, 5)) +
  labs(x = NULL, y = NULL,
       title = "Genome \u00d7 class heatmap (124 genomes \u00d7 21 classes)") +
  theme_pub(base_size = 8) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# --- Panel C: Novelty ---------------------------------------------------------
novelty_sum <- novelty %>% count(novelty_category) %>%
  mutate(pct = 100 * n / sum(n))

p3c <- ggplot(novelty_sum, aes(x = reorder(novelty_category, -n), y = n,
                                fill = novelty_category)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(
    "putatively_novel" = "#E58601",
    "related"          = "#B2182B",
    "known"            = "#2166AC"
  )) +
  scale_x_discrete(labels = c("putatively_novel" = "Putatively\nnovel",
                               "related"          = "Related",
                               "known"            = "Known")) +
  labs(x = NULL, y = "Number of BGCs",
       title = "MIBiG-based novelty classification") +
  theme_pub() +
  expand_limits(y = max(novelty_sum$n) * 1.15)

# --- Combine + save ----------------------------------------------------------
fig3 <- (p3a / p3c | p3b) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "Fig3_biosynthetic_diversity.pdf"), fig3,
       width = 14, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "Fig3_biosynthetic_diversity.png"), fig3,
       width = 14, height = 9, dpi = 300)

cat("Fig 3 saved:", fig_dir, "\n")
