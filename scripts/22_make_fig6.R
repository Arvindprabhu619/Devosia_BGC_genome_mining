#!/usr/bin/env Rscript
# =============================================================================
# 22_make_fig6.R
# Figure 6 — Ecological + Evolutionary Structure
#   Panel A: Niche × BGC class association heatmap (log2 OR + FDR stars)
#   Panel B: GCF evolutionary age distribution
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

# --- Paths -------------------------------------------------------------------
base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
in_dir   <- file.path(base_dir, "results/11_stats")
out_dir  <- file.path(base_dir, "figures/Figure6")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read data ---------------------------------------------------------------
assoc <- read.delim(file.path(in_dir, "ecological_association_classes.tsv"),
                    stringsAsFactors = FALSE)
ages  <- read.delim(file.path(in_dir, "gcf_ages.tsv"),
                    stringsAsFactors = FALSE)

# --- Panel A: Niche × class heatmap (log2 OR + significance stars) -----------
assoc_a <- assoc %>%
  mutate(
    log2_or = log2(pmax(odds_ratio, 0.01)),
    sig     = ifelse(fdr < 0.001, "***",
              ifelse(fdr < 0.01,  "**",
              ifelse(fdr < 0.05,  "*", "")))
  )

# Order classes by number of significant niches (most significant on top)
class_order <- assoc_a %>%
  group_by(class) %>%
  summarise(n_sig = sum(fdr < 0.05), .groups = "drop") %>%
  arrange(desc(n_sig)) %>%
  pull(class)
assoc_a$class <- factor(assoc_a$class, levels = rev(class_order))

# Order niches alphabetically
niche_order <- sort(unique(assoc_a$niche))
assoc_a$niche <- factor(assoc_a$niche, levels = niche_order)

pA <- ggplot(assoc_a, aes(x = niche, y = class, fill = log2_or)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low      = "#2166AC",       # blue  = depleted (OR < 1)
    mid      = "white",          # white = neutral (OR ≈ 1)
    high     = "#B2182B",        # red   = enriched (OR > 1)
    midpoint = 0,
    limits   = c(-3, 3),
    oob      = scales::squish,
    name     = "log2(OR)"
  ) +
  labs(
    title = "A  Niche × BGC class associations",
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid        = element_blank(),
    axis.text.x       = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y       = element_text(size = 9),
    plot.title        = element_text(face = "bold", size = 12),
    legend.position   = "right",
    legend.key.height = unit(1.2, "cm")
  )

# --- Panel B: GCF evolutionary ages -----------------------------------------
age_levels <- c("ancient", "old", "intermediate", "recent", "singleton")
ages$age_category <- factor(ages$age_category, levels = age_levels)

age_summary <- ages %>%
  count(age_category, .drop = FALSE)

age_colors <- c(
  ancient      = "#1B7837",
  old          = "#5AAE61",
  intermediate = "#A6DBA0",
  recent       = "#E7D4E8",
  singleton    = "#9970AB"
)

pB <- ggplot(age_summary, aes(x = age_category, y = n, fill = age_category)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = age_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "B  GCF evolutionary ages",
    x     = NULL,
    y     = "Number of GCFs"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 11),
    axis.text.y        = element_text(size = 10),
    plot.title         = element_text(face = "bold", size = 12)
  )

# --- Combine -----------------------------------------------------------------
fig6 <- pA / pB + plot_layout(heights = c(2, 1))

# --- Save --------------------------------------------------------------------
ggsave(file.path(out_dir, "Figure6.pdf"), fig6,
       width = 10, height = 9)
ggsave(file.path(out_dir, "Figure6.png"), fig6,
       width = 10, height = 9, dpi = 300)

cat("Figure 6 saved to:", out_dir, "\n")
cat("  - Figure6.pdf\n")
cat("  - Figure6.png\n")
