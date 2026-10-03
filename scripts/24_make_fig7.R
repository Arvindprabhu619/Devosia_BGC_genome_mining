#!/usr/bin/env Rscript
# =============================================================================
# 24_make_fig7.R
# Figure 7 — Prioritized BGC candidates
#   Panel A: Product class composition of top 7 candidates
#   Panel B: Score breakdown per candidate
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

# --- Paths -------------------------------------------------------------------
base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
in_file  <- file.path(base_dir, "results/12_prioritization/prioritized_candidates.tsv")
out_dir  <- file.path(base_dir, "figures/Figure7")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read data ---------------------------------------------------------------
cand <- read.delim(in_file, stringsAsFactors = FALSE)

# Parse product classes (split on ;)
classes_long <- cand %>%
  select(rank, bgc_id, genome_accession, product_classes) %>%
  separate_rows(product_classes, sep = ";") %>%
  mutate(
    product_classes = trimws(product_classes),
    candidate_label = paste0(rank, ". ", genome_accession)
  )

# Parse reasons into scoring components
reasons_long <- cand %>%
  select(rank, genome_accession, reasons) %>%
  separate_rows(reasons, sep = ";") %>%
  mutate(
    component = sub("\\+.*", "", reasons),
    value     = as.numeric(sub(".*\\+", "", reasons)),
    candidate_label = paste0(rank, ". ", genome_accession)
  ) %>%
  filter(!is.na(value))

# --- Panel A: product class composition --------------------------------------
# Order candidates by rank (1 at top)
candidate_order <- cand$rank
label_order <- paste0(cand$rank, ". ", cand$genome_accession)
classes_long$candidate_label <- factor(classes_long$candidate_label,
                                       levels = rev(label_order))

# Color palette (colorblind-safe, ~15 classes)
class_colors <- c(
  "RiPP-like"          = "#1F77B4",
  "hydrogen-cyanide"   = "#AEC7E8",
  "NRPS-like"          = "#FF7F0E",
  "T1PKS"              = "#FFBB78",
  "T3PKS"              = "#2CA02C",
  "hserlactone"        = "#98DF8A",
  "lassopeptide"       = "#D62728",
  "RRE-containing"     = "#FF9896"
)

pA <- ggplot(classes_long,
             aes(x = candidate_label, y = product_classes,
                 fill = product_classes)) +
  geom_tile(color = "white", linewidth = 0.5, width = 0.9, height = 0.7) +
  geom_text(aes(label = product_classes), size = 3, color = "white",
            fontface = "bold") +
  scale_fill_manual(values = class_colors, guide = "none") +
  labs(
    title = "A  Product classes of top 7 prioritized candidates",
    x     = NULL,
    y     = NULL
  ) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid        = element_blank(),
    axis.text.x       = element_text(size = 9),
    axis.text.y       = element_text(size = 9),
    plot.title        = element_text(face = "bold", size = 12)
  )

# --- Panel B: score breakdown ------------------------------------------------
# Order components
comp_levels <- c("novelty", "hybrid", "class", "candidates")
reasons_long$component <- factor(reasons_long$component, levels = comp_levels)
reasons_long$candidate_label <- factor(reasons_long$candidate_label,
                                        levels = rev(label_order))

comp_colors <- c(
  novelty    = "#B2182B",
  hybrid     = "#EF8A62",
  class      = "#67A9CF",
  candidates = "#2166AC"
)

pB <- ggplot(reasons_long,
             aes(x = candidate_label, y = value, fill = component)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = comp_colors,
                    name = "Scoring component",
                    breaks = comp_levels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "B  Score breakdown per candidate",
    x     = NULL,
    y     = "Score"
  ) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x        = element_text(size = 9),
    axis.text.y        = element_text(size = 9),
    plot.title         = element_text(face = "bold", size = 12),
    legend.position    = "bottom"
  )

# --- Combine -----------------------------------------------------------------
fig7 <- pA / pB +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 7  Prioritized BGC candidates",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# --- Save --------------------------------------------------------------------
ggsave(file.path(out_dir, "Figure7.pdf"), fig7,
       width = 11, height = 8)
ggsave(file.path(out_dir, "Figure7.png"), fig7,
       width = 11, height = 8, dpi = 300)

cat("Figure 7 saved to:", out_dir, "\n")
cat("  - Figure7.pdf\n")
cat("  - Figure7.png\n")
