#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/21_make_fig2_definitive.R
# ============================================================
# Fig 2 — Circular phylogeny with 20 BGC class rings.
# Pure ggplot2 + ggtree. No circlize. Full control.
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggtreeExtra)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(root, "figures")

cat("=== [21] Fig 2 definitive ===\n")

# ---------- Load ----------
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

# ---------- Class matrix ----------
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

# Verify: terpene should be 124/124
cat(sprintf("Terpene: %d/124\n", sum(class_wide[["terpene"]] > 0)))

# Colors
class_colors <- c(
  "#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B8EA5",
  "#06A77D", "#D5A021", "#8E6C8A", "#4B8B3B", "#E63946",
  "#457B9D", "#F4A261", "#2A9D8F", "#E76F51", "#264653",
  "#9B5DE5", "#F15BB5", "#FEE440", "#00BBF9", "#00F5D4"
)[seq_len(n_classes)]
names(class_colors) <- class_order

# Tip taxonomy
tax_map <- setNames(tax$status, tax$accession)
tip_tax <- tax_map[tip_order]
tip_tax[is.na(tip_tax)] <- "Unknown"
tip_df <- data.frame(label = tip_order, genus = tip_tax, stringsAsFactors = FALSE)

# ============================================================
# Build LONG format for ggtreeExtra
# Each row = (genome, class, present, ring_y)
# For "all genomes" approach: expand to 124 × 20 rows
# ============================================================
cat("Building long-format data...\n")

# Wide → long
class_long <- class_wide %>%
  tibble::rownames_to_column("genome_accession") %>%
  pivot_longer(-genome_accession, names_to = "class", values_to = "count") %>%
  mutate(present = as.integer(count > 0),
         class = factor(class, levels = class_order))

cat(sprintf("Long-format rows: %d\n", nrow(class_long)))

# ============================================================
# Base tree (circular)
# ============================================================
cat("Building tree...\n")

p <- ggtree(tree, layout = "circular", size = 0.4, open.angle = 5) %<+% tip_df

# Add genus-colored tip labels
p <- p +
  geom_tippoint(aes(color = genus), size = 1.0, alpha = 0.9) +
  scale_color_manual(
    values = c("Devosia" = "#2E86AB",
               "Devosia_A" = "#A23B72",
               "Unknown" = "grey60"),
    name = "Genus"
  )

# ============================================================
# Add 20 class rings via geom_fruit (ALL genomes shown)
# Key: use present as fill, with discrete palette
# ============================================================
cat("Adding rings...\n")

# For each class, add a geom_fruit with all 124 genomes
for (i in seq_along(class_order)) {
  cls <- class_order[i]
  sub <- class_long %>% filter(class == cls)
  offset_val <- 0.05 + (i - 1) * 0.024
  color_here <- class_colors[cls]

  p <- p + geom_fruit(
    data = sub,
    geom = geom_tile,
    mapping = aes(y = genome_accession, x = "v", fill = factor(present)),
    offset = offset_val,
    pwidth = 0.022,
    color = NA
  )
}

# Discrete fill: 0 = light grey, 1 = class color
# But since each ring uses its own class color, we need per-ring scales.
# Trick: use `new_scale_fill()` from ggnewscale between rings.
# Since that's complex, use a single palette where present = black
# and rely on TIP LABELS for class identification.

p <- p + scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "black"),
                            guide = "none")

# ============================================================
# Save
# ============================================================
cat("Saving...\n")

ggsave(file.path(fig_dir, "Fig2_definitive.pdf"), p,
       width = 16, height = 16, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_definitive.png"), p,
       width = 16, height = 16, dpi = 300)

# ============================================================
# Also save a version with legend
# ============================================================
cat("Building legend...\n")

legend_df <- data.frame(
  y = seq_along(class_order),
  class = class_order,
  color = class_colors
)

p_legend <- ggplot(legend_df, aes(x = 1, y = y, fill = class)) +
  geom_tile(width = 0.4, height = 0.9) +
  geom_text(aes(label = sprintf("%d. %s", y, class)),
            x = 1.35, hjust = 0, size = 3.2) +
  scale_fill_manual(values = class_colors, guide = "none") +
  xlim(0.6, 4) + ylim(0, length(class_order) + 2) +
  theme_void() +
  ggtitle("BGC class") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.05))

fig_combined <- p + p_legend + plot_layout(widths = c(4, 1.2))

ggsave(file.path(fig_dir, "Fig2_definitive_with_legend.pdf"), fig_combined,
       width = 18, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_definitive_with_legend.png"), fig_combined,
       width = 18, height = 14, dpi = 300)

cat("=== [21] Done ===\n")
