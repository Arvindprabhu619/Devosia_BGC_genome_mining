#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/23_fig2_streptomyces_style.R
# ============================================================
# Fig 2 — Streptomyces-style circular phylogeny with 20 class rings.
# Uses ggtree + ggtreeExtra + ggnewscale for correct multi-color rings.
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(root, "figures")

cat("=== [23] Fig 2 Streptomyces style ===\n")

# ---------- Load data ----------
tree <- read.tree(file.path(root, "results/05b_phylogeny/devosia_pruned.tree"))
classes <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv"),
                       stringsAsFactors = FALSE)
tax <- read.delim(file.path(root, "results/05_gtdbtk/taxonomy_assignments.tsv"),
                   stringsAsFactors = FALSE)

if (!is.binary(tree)) tree <- multi2di(tree)
tree$edge.length[tree$edge.length == 0] <- 1e-6
tip_order <- tree$tip.label
n_tips <- length(tip_order)
cat(sprintf("Tips: %d\n", n_tips))

# ---------- Class presence matrix ----------
class_wide <- classes %>%
  count(genome_accession, product_class, name = "n") %>%
  pivot_wider(id_cols = genome_accession, names_from = product_class,
              values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(class_wide) <- class_wide$genome_accession
class_wide$genome_accession <- NULL
class_wide <- class_wide[tip_order, , drop = FALSE]

class_order <- names(sort(colSums(class_wide > 0), decreasing = TRUE))
class_wide <- class_wide[, class_order, drop = FALSE]
n_classes <- length(class_order)
cat(sprintf("Classes: %d\n", n_classes))

# ---------- Class colors ----------
# Distinct, publication-quality palette for 20 classes
class_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
  "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5"
)[seq_len(n_classes)]
names(class_colors) <- class_order

# ---------- Long-format: one row per (genome, class) ----------
class_long <- class_wide %>%
  tibble::rownames_to_column("genome") %>%
  pivot_longer(-genome, names_to = "class", values_to = "count") %>%
  mutate(present = as.integer(count > 0),
         class = factor(class, levels = class_order))

# ---------- Tip taxonomy for tip colors ----------
tax_map <- setNames(tax$status, tax$accession)
tip_tax <- tax_map[tip_order]
tip_tax[is.na(tip_tax)] <- "Unknown"
tip_df <- data.frame(label = tip_order, genus = tip_tax, stringsAsFactors = FALSE)

# ============================================================
# Build ggtree base
# ============================================================
cat("Building tree base...\n")

p <- ggtree(tree, layout = "circular", size = 0.35, open.angle = 5) %<+% tip_df

# Color tips by genus
p <- p + geom_tippoint(aes(color = genus), size = 0.8, alpha = 0.9) +
  scale_color_manual(
    values = c("Devosia" = "#1F77B4", "Devosia_A" = "#D62728", "Unknown" = "grey60"),
    name = "Genus"
  )

# ============================================================
# Add 20 class rings — each with its own color via new_scale_fill
# ============================================================
cat("Adding", n_classes, "rings with new_scale_fill between each...\n")

for (i in seq_len(n_classes)) {
  cls <- class_order[i]
  sub <- class_long %>% filter(class == cls)
  color_i <- class_colors[i]

  # Add ring — present = class color, absent = light grey
  p <- p + geom_fruit(
    data = sub,
    geom = geom_tile,
    mapping = aes(y = genome, x = "v", fill = factor(present)),
    offset = 0.03 + (i - 1) * 0.022,
    pwidth = 0.018,
    color = NA
  ) +
    scale_fill_manual(
      values = c("0" = "#F0F0F0", "1" = color_i),
      guide = "none"
    ) +
    ggnewscale::new_scale_fill()   # ← critical: reset fill scale for next ring
}

# ============================================================
# Save plot (tree only, no legend)
# ============================================================
cat("Saving tree + rings...\n")
ggsave(file.path(fig_dir, "Fig2_streptomyces_style.pdf"), p,
       width = 14, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_streptomyces_style.png"), p,
       width = 14, height = 14, dpi = 300)

# ============================================================
# Build clean legend panel
# ============================================================
cat("Building legend panel...\n")

legend_df <- data.frame(
  class = factor(class_order, levels = rev(class_order)),
  y = seq_along(class_order),
  color = rev(class_colors[class_order]),
  label = sprintf("%d. %s", rev(seq_along(class_order)), rev(class_order))
)

p_legend <- ggplot(legend_df, aes(x = 1, y = y, fill = class)) +
  geom_tile(width = 0.35, height = 0.85) +
  geom_text(aes(label = label), x = 1.3, hjust = 0, size = 3.5) +
  scale_fill_manual(values = rev(class_colors), guide = "none") +
  xlim(0.7, 4.5) +
  ylim(0.5, n_classes + 0.5) +
  theme_void() +
  ggtitle("BGC class") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.02))

# Combine tree + legend side-by-side
fig_combined <- p + p_legend +
  plot_layout(widths = c(4, 1.2)) +
  plot_annotation(
    title = "Figure 2. Phylogenetic distribution of 20 BGC classes across 124 Devosia genomes",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0)
    )
  )

ggsave(file.path(fig_dir, "Fig2_with_legend.pdf"), fig_combined,
       width = 18, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_with_legend.png"), fig_combined,
       width = 18, height = 14, dpi = 300)

cat("=== [23] Done ===\n")
cat("Saved:\n")
cat("  figures/Fig2_streptomyces_style.pdf/png\n")
cat("  figures/Fig2_with_legend.pdf/png\n")
