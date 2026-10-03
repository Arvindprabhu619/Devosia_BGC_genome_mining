#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/17_make_class_heatmap_circlize.R
# ============================================================
# Purpose: Dense circular heatmap showing BGC class presence
#          across the 124-genome phylogeny.
#
# Style: matches Streptomyces-style concentric ring figures.
#
# Outputs:
#   figures/Fig2_class_heatmap_ring.pdf + .png
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))

tree_file <- file.path(root, "results/05b_phylogeny/devosia_pruned.tree")
classes_file <- file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv")
novelty_file <- file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv")
tax_file <- file.path(root, "results/05_gtdbtk/taxonomy_assignments.tsv")
fig_dir <- file.path(root, "figures")

cat("=== [17] Circlize class heatmap ring ===\n")

# ---------- Load tree ----------
tree <- read.tree(tree_file)
tree <- multi2di(tree)
tree$edge.length[tree$edge.length == 0] <- 1e-6
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

# ---------- Load classes ----------
classes <- read.delim(classes_file, stringsAsFactors = FALSE)

# Class × genome matrix
class_counts <- classes %>%
  count(genome_accession, product_class, name = "n") %>%
  complete(genome_accession = tree$tip.label,
           product_class = unique(classes$product_class),
           fill = list(n = 0)) %>%
  pivot_wider(id_cols = genome_accession, names_from = product_class,
              values_from = n, values_fill = 0) %>%
  as.data.frame()

# Reorder to match tree
rownames(class_counts) <- class_counts$genome_accession
class_counts$genome_accession <- NULL
class_counts <- class_counts[tree$tip.label, , drop = FALSE]

# Order classes by total presence
class_totals <- colSums(class_counts > 0)
class_order <- names(sort(class_totals, decreasing = TRUE))
class_counts <- class_counts[, class_order, drop = FALSE]

cat(sprintf("Classes: %d\n", ncol(class_counts)))

# ---------- Load novelty ----------
novelty <- read.delim(novelty_file, stringsAsFactors = FALSE)
gnovel <- novelty %>%
  group_by(genome_accession) %>%
  summarise(prop_novel = sum(novelty_category == "putatively_novel") / n(),
            .groups = "drop")
gnovel <- gnovel[match(tree$tip.label, gnovel$genome_accession), ]

# ---------- Load taxonomy ----------
tax <- read.delim(tax_file, stringsAsFactors = FALSE)
tax_map <- setNames(tax$status, tax$accession)
tax_status <- tax_map[tree$tip.label]
tax_status[is.na(tax_status)] <- "Unknown"

# ---------- Color palettes ----------
# 21 class colors - use Set3 + Spectral blend
n_classes <- ncol(class_counts)
class_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_classes)
names(class_colors) <- class_order

# Tax palette
tax_colors <- c("Devosia" = "#2E86AB",
                "Devosia_A" = "#A23B72",
                "Unknown" = "grey70")

# ---------- Open PDF device ----------
pdf(file.path(fig_dir, "Fig2_class_heatmap_ring.pdf"),
    width = 14, height = 14)
png(file.path(fig_dir, "Fig2_class_heatmap_ring.png"),
    width = 14, height = 14, units = "in", res = 300)

par(mar = c(1, 1, 1, 1))

# ---------- Circular layout ----------
circos.par(
  start.degree = 90,
  gap.degree = 0.5,
  track.margin = c(0.005, 0.005),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = "tree", xlim = c(0, length(tree$tip.label)))

# ---------- Add tree as first track ----------
circos.track(ylim = c(0, 1), track.height = 0.15,
             bg.border = NA, panel.fun = function(x, y) {
  # Place tip labels around the circle
  circos.text(CELL_META$xcenter, 0.5, "", facing = "clockwise",
              niceFacing = TRUE, cex = 0.5)
})

# ---------- Class presence rings ----------
# Each class gets its own narrow ring
for (i in seq_len(n_classes)) {
  cls <- class_order[i]
  values <- class_counts[[cls]]
  presence <- as.integer(values > 0)

  circos.track(
    ylim = c(0, 1),
    track.height = 0.02,
    bg.border = "grey90",
    bg.col = "white",
    panel.fun = function(x, y) {
      for (j in seq_along(presence)) {
        if (presence[j] == 1) {
          x_start <- j - 1
          x_end <- j
          circos.rect(x_start, 0, x_end, 1,
                      col = class_colors[cls],
                      border = NA)
        }
      }
    }
  )
}

# ---------- Add legend tracks on outer edge ----------
# Genome abundance ring (outermost)
track_abundance <- rowSums(class_counts)
max_abund <- max(track_abundance)

circos.track(
  ylim = c(0, max_abund),
  track.height = 0.05,
  bg.border = NA,
  bg.col = NA,
  panel.fun = function(x, y) {
    for (j in seq_along(track_abundance)) {
      circos.rect(j - 1, 0, j, track_abundance[j],
                  col = "#2E86AB", border = NA)
    }
  }
)

# Novelty proportion ring (outermost)
circos.track(
  ylim = c(0, 1),
  track.height = 0.03,
  bg.border = "grey90",
  panel.fun = function(x, y) {
    for (j in seq_len(nrow(gnovel))) {
      prop <- gnovel$prop_novel[j]
      if (is.na(prop)) prop <- 0
      col <- colorRampPalette(c("white", "#C73E1D"))(100)[max(1, round(prop * 99) + 1)]
      circos.rect(j - 1, 0, j, 1, col = col, border = NA)
    }
  }
)

# ---------- Tree as outermost arc labels ----------
# Draw the tree structure as a set of segments on the outside
circos.track(
  ylim = c(0, 1),
  track.height = 0.06,
  bg.border = NA,
  panel.fun = function(x, y) {
    for (j in seq_along(tax_status)) {
      col <- tax_colors[tax_status[j]]
      if (is.na(col)) col <- "grey70"
      circos.rect(j - 1, 0, j, 1, col = col, border = "white", lwd = 0.3)
    }
  }

)
# Close devices
circos.clear()
dev.off()  # close PDF
dev.off()  # close PNG

cat("  Saved Fig2_class_heatmap_ring.pdf + .png\n")
cat("=== [17] Done ===\n")
