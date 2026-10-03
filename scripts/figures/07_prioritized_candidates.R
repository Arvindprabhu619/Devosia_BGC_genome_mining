#!/usr/bin/env Rscript
# =============================================================================
# 24_make_fig7.R
# Figure 7 — Prioritized BGC candidates
#   Panel A: Product class composition of top 7 candidates (split-row per candidate)
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

# Build candidate label
cand$candidate_label <- paste0(cand$rank, ". ", cand$genome_accession)
label_order <- cand$candidate_label

# --- Parse product classes ---------------------------------------------------
# Each candidate has 2 classes separated by ";"
classes_long <- cand %>%
  select(rank, candidate_label, product_classes) %>%
  separate_rows(product_classes, sep = ";") %>%
  mutate(
    product_classes = trimws(product_classes),
    # Position within candidate (1 or 2) for split-row layout
    position = ave(seq_along(product_classes), candidate_label, FUN = seq_along)
  )

# Order candidates by rank
classes_long$candidate_label <- factor(classes_long$candidate_label,
                                       levels = rev(label_order))

# Color palette (colorblind-safe)
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

# --- Panel A: split-row per candidate ----------------------------------------
# x = position (1 or 2), y = candidate
pA <- ggplot(classes_long,
             aes(x = factor(position), y = candidate_label,
                 fill = product_classes)) +
  geom_tile(color = "white", linewidth = 1, width = 0.95, height = 0.85) +
  geom_text(aes(label = product_classes),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = class_colors, guide = "none") +
  scale_x_discrete(
    labels = c("Class 1", "Class 2"),
    expand = expansion(add = 0.3)
  ) +
  labs(
    title = "A  Product classes of top 7 prioritized candidates",
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid        = element_blank(),
    axis.text.x       = element_text(size = 10),
    axis.text.y       = element_text(size = 10),
    plot.title        = element_text(face = "bold", size = 12),
    plot.margin       = margin(5, 25, 5, 5)
  )

# --- Panel B: score breakdown ------------------------------------------------
# Parse reasons: "novelty+3;hybrid+1;class+1;candidates+1"
reasons_long <- cand %>%
  select(candidate_label, reasons) %>%
  separate_rows(reasons, sep = ";") %>%
  mutate(
    component = sub("\\+.*", "", reasons),
    value     = as.numeric(sub(".*\\+", "", reasons))
  ) %>%
  filter(!is.na(value), value > 0)

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
  geom_col(width = 0.7, color = "white", linewidth = 0.4) +
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5),
            size = 3.8, color = "white", fontface = "bold") +
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
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 10),
    axis.text.y        = element_text(size = 10),
    plot.title         = element_text(face = "bold", size = 12),
    legend.position    = "bottom",
    legend.title       = element_text(size = 10, face = "bold"),
    legend.text        = element_text(size = 9),
    plot.margin        = margin(5, 5, 5, 5)
  )

# --- Combine -----------------------------------------------------------------
fig7 <- pA / pB +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 7  Prioritized BGC candidates",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0)
    )
  )

# --- Save --------------------------------------------------------------------
ggsave(file.path(out_dir, "Figure7.pdf"), fig7,
       width = 11, height = 9)
ggsave(file.path(out_dir, "Figure7.png"), fig7,
       width = 11, height = 9, dpi = 300)

cat("Figure 7 saved to:", out_dir, "\n")
cat("  - Figure7.pdf\n")
cat("  - Figure7.png\n")
