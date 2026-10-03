#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/20_make_fig2_final.R
# ============================================================
# Publication-grade Fig 2 — Streptomyces-style circular phylogeny.
# ============================================================

suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggtreeExtra)
  library(ggplot2); library(dplyr); library(tidyr)
  library(RColorBrewer); library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(root, "figures")

cat("=== [20] Fig 2 final ===\n")

# ---------- Load data ----------
tree <- read.tree(file.path(root, "results/05b_phylogeny/devosia_pruned.tree"))
classes <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv"),
                       stringsAsFactors = FALSE)
tax <- read.delim(file.path(root, "results/05_gtdbtk/taxonomy_assignments.tsv"),
                   stringsAsFactors = FALSE)

cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

# ---------- Class presence matrix ----------
class_wide <- classes %>%
  count(genome_accession, product_class, name = "n") %>%
  complete(genome_accession = tree$tip.label,
           product_class = unique(classes$product_class),
           fill = list(n = 0)) %>%
  pivot_wider(id_cols = genome_accession,
              names_from = product_class, values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(class_wide) <- class_wide$genome_accession
class_wide$genome_accession <- NULL
class_wide <- class_wide[tree$tip.label, , drop = FALSE]

class_order <- names(sort(colSums(class_wide > 0), decreasing = TRUE))
class_wide <- class_wide[, class_order, drop = FALSE]
n_classes <- length(class_order)
cat(sprintf("BGC classes: %d\n", n_classes))

# ---------- Class colors ----------
class_colors <- c(
  "#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B8EA5",
  "#06A77D", "#D5A021", "#8E6C8A", "#4B8B3B", "#E63946",
  "#457B9D", "#F4A261", "#2A9D8F", "#E76F51", "#264653",
  "#9B5DE5", "#F15BB5", "#FEE440", "#00BBF9", "#00F5D4",
  "#8E7DBE"
)[seq_len(n_classes)]
names(class_colors) <- class_order

# ---------- Tip taxonomy ----------
tax_map <- setNames(tax$status, tax$accession)
tip_tax <- tax_map[tree$tip.label]
tip_tax[is.na(tip_tax)] <- "Unknown"
tip_df <- data.frame(label = tree$tip.label, genus = tip_tax, stringsAsFactors = FALSE)

# ---------- Long format ----------
class_long <- class_wide %>%
  tibble::rownames_to_column("genome_accession") %>%
  pivot_longer(-genome_accession, names_to = "class", values_to = "count") %>%
  mutate(present = as.integer(count > 0),
         class = factor(class, levels = class_order))

# ---------- Build ggtree ----------
cat("Building tree...\n")

p <- ggtree(tree, layout = "circular", size = 0.3, open.angle = 15) %<+% tip_df

# Add 21 rings with distinct colors
cat("Adding 21 rings...\n")
for (i in seq_along(class_order)) {
  cls <- class_order[i]
  sub <- class_long %>% filter(class == cls) %>% mutate(ring_id = cls)
  p <- p + geom_fruit(
    data = sub,
    geom = geom_tile,
    mapping = aes(y = genome_accession, fill = ring_id),
    offset = 0.02 + (i - 1) * 0.022,
    pwidth = 0.018,
    color = NA,
    show.legend = FALSE
  )
}

# Fill scale for class colors
p <- p + scale_fill_manual(values = class_colors, guide = "none")

# ---------- Save ----------
cat("Saving...\n")
ggsave(file.path(fig_dir, "Fig2_final.pdf"), p,
       width = 14, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_final.png"), p,
       width = 14, height = 14, dpi = 300)

cat("=== [20] Done ===\n")
