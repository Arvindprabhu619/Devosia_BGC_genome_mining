#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/21_make_fig2_definitive.R
# ============================================================
# Fig 2 — Streptomyces-style circular phylogeny.
# Includes: real tree, 20 class rings, scale bar, panel title.
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree); library(ggtreeExtra)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(root, "figures")

cat("=== [21] Fig 2 polished ===\n")

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
  complete(genome_accession = tip_order,
           product_class = unique(classes$product_class),
           fill = list(n = 0)) %>%
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

# Colors
class_colors <- c(
  "#2E86AB","#A23B72","#F18F01","#C73E1D","#3B8EA5",
  "#06A77D","#D5A021","#8E6C8A","#4B8B3B","#E63946",
  "#457B9D","#F4A261","#2A9D8F","#E76F51","#264653",
  "#9B5DE5","#F15BB5","#FEE440","#00BBF9","#00F5D4"
)[seq_len(n_classes)]
names(class_colors) <- class_order

# Tip taxonomy
tax_map <- setNames(tax$status, tax$accession)
tip_tax <- tax_map[tip_order]
tip_tax[is.na(tip_tax)] <- "Unknown"

tip_df <- data.frame(label = tip_order, genus = tip_tax, stringsAsFactors = FALSE)

# ---------- Build tree base with ggtree ----------
cat("Building tree...\n")

p_tree <- ggtree(tree, layout = "circular", size = 0.35, open.angle = 5) %<+% tip_df

# Color tip labels by genus
p_tree <- p_tree +
  geom_tiplab(aes(color = genus), size = 1.5, offset = 0.5, align = FALSE) +
  scale_color_manual(values = c("Devosia" = "#2E86AB",
                                 "Devosia_A" = "#A23B72",
                                 "Unknown" = "grey70"),
                     name = "Genus")

# ---------- Long format for rings ----------
class_long <- class_wide %>%
  tibble::rownames_to_column("genome_accession") %>%
  pivot_longer(-genome_accession, names_to = "class", values_to = "count") %>%
  mutate(present = as.integer(count > 0),
         class = factor(class, levels = class_order))

# ---------- Add 20 rings ----------
cat("Adding 20 rings...\n")
ring_offset_base <- 0.10
ring_spacing <- 0.022
ring_width <- 0.018

for (i in seq_along(class_order)) {
  cls <- class_order[i]
  sub <- class_long %>% filter(class == cls)

  p_tree <- p_tree + geom_fruit(
    data = sub,
    geom = geom_tile,
    mapping = aes(y = genome_accession, fill = factor(present)),
    offset = ring_offset_base + (i - 1) * ring_spacing,
    pwidth = ring_width,
    color = NA
  )
}

# Discrete fill: 0 = light grey, 1 = black
p_tree <- p_tree + scale_fill_manual(
  values = c("0" = "#F0F0F0", "1" = "black"),
  guide = "none"
)

# ---------- Save PDF + PNG ----------
cat("Saving...\n")
ggsave(file.path(fig_dir, "Fig2_definitive.pdf"), p_tree,
       width = 16, height = 16, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_definitive.png"), p_tree,
       width = 16, height = 16, dpi = 300)

# ============================================================
# Also save a version with legend via patchwork
# ============================================================
cat("Building legend version...\n")

# Create separate legend plot
legend_df <- data.frame(
  class = factor(class_order, levels = rev(class_order)),
  y = seq_along(class_order),
  color = class_colors
)

p_legend <- ggplot(legend_df, aes(x = 1, y = y, fill = class)) +
  geom_tile(width = 0.5, height = 0.85) +
  geom_text(aes(label = class), x = 1.35, hjust = 0, size = 3.5) +
  scale_fill_manual(values = class_colors, guide = "none") +
  xlim(0.7, 3) +
  theme_void() +
  ggtitle("BGC class") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.05),
        plot.margin = margin(10, 10, 10, 10))

fig2_combined <- p_tree | p_legend
fig2_combined <- fig2_combined + plot_layout(widths = c(4, 1))

ggsave(file.path(fig_dir, "Fig2_definitive_with_legend.pdf"), fig2_combined,
       width = 20, height = 16, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_definitive_with_legend.png"), fig2_combined,
       width = 20, height = 16, dpi = 300)

cat("=== [21] Done ===\n")
cat("Saved:\n")
cat("  figures/Fig2_definitive.pdf/png\n")
cat("  figures/Fig2_definitive_with_legend.pdf/png\n")
