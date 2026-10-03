#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/16_make_ring_figure.R
# ============================================================
# Purpose: Circular phylogenetic ring figures showing BGC
#          class distribution and multi-track genome properties.
#
# Outputs:
#   figures/Fig2_class_tree.pdf + .png         (21 class rings)
#   figures/Fig9_multitrack_tree.pdf + .png    (abundance + novelty + GCFs)
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggtreeExtra)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(RColorBrewer)
  library(patchwork)
  library(cowplot)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))

tree_file <- file.path(root, "results/05b_phylogeny/devosia_pruned.tree")
regions_file <- file.path(root, "results/06_antismash/parsed/bgc_regions.tsv")
classes_file <- file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv")
novelty_file <- file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv")
gcf_file <- file.path(root, "results/09_bigscape/parsed/gcf_prevalence_c0.7.tsv")
genome_list <- file.path(root, "results/05_gtdbtk/final_genome_list.txt")
tax_file <- file.path(root, "results/05_gtdbtk/taxonomy_assignments.tsv")

fig_dir <- file.path(root, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== [16] Ring figures ===\n")
cat("Loading data...\n")

# ---------- Load tree ----------
tree <- read.tree(tree_file)
tree <- multi2di(tree)
tree$edge.length[tree$edge.length == 0] <- 1e-6
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

# ---------- Load genomes ----------
genomes <- readLines(genome_list)
genomes <- genomes[genomes != ""]
common <- intersect(tree$tip.label, genomes)
tree <- keep.tip(tree, common)
cat(sprintf("Common tips: %d\n", length(common)))

# ---------- Load taxonomy ----------
tax <- read.delim(tax_file, stringsAsFactors = FALSE)
tax$genome <- tax$accession
tax_map <- setNames(tax$status, tax$genome)

# ---------- Load BGC classes ----------
regions <- read.delim(regions_file, stringsAsFactors = FALSE)
classes <- read.delim(classes_file, stringsAsFactors = FALSE)

# Class-by-genome presence matrix (wide)
class_wide <- classes %>%
  count(genome_accession, product_class, name = "n") %>%
  complete(genome_accession = common, product_class = unique(classes$product_class),
           fill = list(n = 0)) %>%
  pivot_wider(id_cols = genome_accession, names_from = product_class,
              values_from = n, values_fill = 0)

# Reorder by tree tips
class_wide <- class_wide[match(tree$tip.label, class_wide$genome_accession), ]

# Get class order by total presence
class_order <- classes %>% count(product_class, sort = TRUE) %>% pull(product_class)

# ---------- Load novelty ----------
novelty <- read.delim(novelty_file, stringsAsFactors = FALSE)
genome_novelty <- novelty %>%
  group_by(genome_accession) %>%
  summarise(
    n_total = n(),
    n_novel = sum(novelty_category == "putatively_novel"),
    n_related = sum(novelty_category == "related"),
    n_known = sum(novelty_category == "known"),
    prop_novel = n_novel / n_total,
    .groups = "drop"
  ) %>%
  right_join(data.frame(genome_accession = common), by = "genome_accession") %>%
  replace_na(list(n_total = 0, n_novel = 0, prop_novel = 0))
genome_novelty <- genome_novelty[match(tree$tip.label, genome_novelty$genome_accession), ]

# ---------- Load GCF prevalence ----------
gcf <- read.delim(gcf_file, stringsAsFactors = FALSE)
top_gcfs <- gcf %>% arrange(desc(prevalence)) %>% head(20) %>% pull(gcf_id)

# Build GCF × genome matrix for top 20 GCFs
gcf_wide <- data.frame(genome_accession = tree$tip.label, stringsAsFactors = FALSE)
for (g in top_gcfs) {
  row <- gcf %>% filter(gcf_id == g)
  genomes_with <- strsplit(row$genomes[1], ",")[[1]]
  gcf_wide[[g]] <- as.integer(tree$tip.label %in% genomes_with)
}

# ---------- Palette ----------
# 21 BGC class colors (spectrum)
n_classes <- length(unique(classes$product_class))
class_pal <- colorRampPalette(brewer.pal(11, "Spectral"))(n_classes)
names(class_pal) <- unique(classes$product_class)

# Devosia vs Devosia_A palette
tax_pal <- c("Devosia" = "#2E86AB", "Devosia_A" = "#A23B72")

# ============================================================
# Figure 2 — BGC CLASS RING TREE (Streptomyces-style)
# ============================================================
cat("Generating Fig 2: Class ring tree...\n")

# Prepare long-format data for geom_fruit_list
class_long <- class_wide %>%
  pivot_longer(cols = -genome_accession, names_to = "class", values_to = "count") %>%
  mutate(present = as.integer(count > 0))

# Root tree for visualization
p_tree <- ggtree(tree, layout = "circular", linewidth = 0.3) %<+%
  data.frame(label = tree$tip.label,
             tax_status = tax_map[tree$tip.label])

# Build class rings
class_rings <- list()
for (i in seq_along(class_order)) {
  cls <- class_order[i]
  sub <- class_long %>% filter(class == cls)
  class_rings[[i]] <- geom_fruit(
    data = sub,
    geom = geom_tile,
    mapping = aes(y = genome_accession, x = "present", fill = factor(present)),
    offset = 0.02 + 0.02 * i,
    pwidth = 0.015,
    show.legend = FALSE
  )
}

# Combine base tree + rings
p2 <- p_tree
for (r in class_rings) {
  p2 <- p2 + r
}

# Custom fill for present/absent
p2 <- p2 + scale_fill_manual(values = c("0" = "grey95", "1" = "black"), guide = "none")

# Save
ggsave(file.path(fig_dir, "Fig2_class_tree.pdf"), p2,
       width = 12, height = 12, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_class_tree.png"), p2,
       width = 12, height = 12, dpi = 300)

cat("  Saved Fig2_class_tree\n")

# ============================================================

# ============================================================
# Figure 9 — MULTI-TRACK RING TREE (2 tracks: abundance + novelty)
# ============================================================
cat("Generating Fig 9: Multi-track tree...\n")

# Tracks
track_abundance <- data.frame(
  genome_accession = tree$tip.label,
  n_bgcs = as.integer(table(regions$genome_accession)[tree$tip.label])
)
track_abundance$n_bgcs[is.na(track_abundance$n_bgcs)] <- 0

track_novelty <- data.frame(
  genome_accession = tree$tip.label,
  prop_novel = genome_novelty$prop_novel
)

track_tax <- data.frame(
  genome_accession = tree$tip.label,
  tax_status = tax_map[tree$tip.label]
)
track_tax$tax_status[is.na(track_tax$tax_status)] <- "Unknown"

# Base tree
p9_base <- ggtree(tree, layout = "circular", linewidth = 0.3, aes(color = tax_status)) %<+% track_tax

# Track 1: BGC abundance bar
p9_base <- p9_base + geom_fruit(
  data = track_abundance,
  geom = geom_col,
  mapping = aes(y = genome_accession, x = n_bgcs),
  fill = "#2E86AB",
  offset = 0.05,
  pwidth = 0.3
)

# Track 2: Novelty proportion heatmap
p9_base <- p9_base + geom_fruit(
  data = track_novelty,
  geom = geom_tile,
  mapping = aes(y = genome_accession, x = "novelty", fill = prop_novel),
  offset = 0.4,
  pwidth = 0.1
)

# Final scales
p9_base <- p9_base +
  scale_color_manual(values = tax_pal, name = "Genus", na.value = "grey70") +
  scale_fill_viridis_c(option = "mako", name = "Novelty\nproportion", direction = -1)

# Save
ggsave(file.path(fig_dir, "Fig9_multitrack_tree.pdf"), p9_base,
       width = 12, height = 12, dpi = 300)
ggsave(file.path(fig_dir, "Fig9_multitrack_tree.png"), p9_base,
       width = 12, height = 12, dpi = 300)

cat("  Saved Fig9_multitrack_tree\n")

cat("=== [16] Done ===\n")
