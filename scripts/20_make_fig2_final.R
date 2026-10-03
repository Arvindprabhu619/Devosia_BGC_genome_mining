#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/20_make_fig2_final.R
# ============================================================
# Publication-grade Fig 2 — Streptomyces-style circular phylogeny
# using circlize for proper presence/absence ring rendering.
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(root, "figures")

cat("=== [20] Fig 2 final (circlize) ===\n")

# ---------- Load data ----------
tree <- read.tree(file.path(root, "results/05b_phylogeny/devosia_pruned.tree"))
classes <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv"),
                       stringsAsFactors = FALSE)
tax <- read.delim(file.path(root, "results/05_gtdbtk/taxonomy_assignments.tsv"),
                   stringsAsFactors = FALSE)

cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

# Make tree binary + add small epsilon to zero-length branches
if (!is.binary(tree)) {
  tree <- multi2di(tree)
  tree$edge.length[tree$edge.length == 0] <- 1e-7
}

# Order tips for plotting
tip_order <- tree$tip.label
n_tips <- length(tip_order)

# ---------- Class presence matrix ----------
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

# Order classes by total presence (most common at innermost ring)
class_order <- names(sort(colSums(class_wide > 0), decreasing = TRUE))
class_wide <- class_wide[, class_order, drop = FALSE]
n_classes <- length(class_order)
cat(sprintf("BGC classes: %d\n", n_classes))

# ---------- Class colors ----------
class_colors <- c(
  "#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B8EA5",
  "#06A77D", "#D5A021", "#8E6C8A", "#4B8B3B", "#E63946",
  "#457B9D", "#F4A261", "#2A9D8F", "#E76F51", "#264653",
  "#9B5DE5", "#F15BB5", "#FEE440", "#00BBF9", "#00F5D4"
)[seq_len(n_classes)]
names(class_colors) <- class_order

# ---------- Tip taxonomy (for tip labels) ----------
tax_map <- setNames(tax$status, tax$accession)
tip_tax <- tax_map[tip_order]
tip_tax[is.na(tip_tax)] <- "Unknown"

# ---------- Init circos ----------
pdf(file.path(fig_dir, "Fig2_final.pdf"), width = 14, height = 14)
png(file.path(fig_dir, "Fig2_final.png"), width = 14, height = 14,
    units = "in", res = 300)

par(mar = c(1, 1, 1, 1))

circos.par(
  start.degree = 90,
  gap.degree = 1,
  track.margin = c(0.002, 0.002),
  points.overflow.warning = FALSE,
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = "genome", xlim = c(0, n_tips))

# ============================================================
# Track 1 (innermost): Phylogenetic tree
# ============================================================
cat("Drawing tree track...\n")

# Compute tree topology as arcs
# We'll draw the tree as a set of straight segments in the ring
# First, get coordinates for each tip
tip_positions <- data.frame(
  tip = tip_order,
  x_start = seq(0, n_tips - 1),
  x_end = seq(1, n_tips),
  stringsAsFactors = FALSE
)

circos.track(
  ylim = c(0, 1),
  track.height = 0.15,
  bg.border = NA,
  panel.fun = function(x, y) {
    # Draw tip labels
    for (i in seq_len(n_tips)) {
      circos.text(i - 0.5, 0.5, tip_order[i],
                  facing = "clockwise",
                  niceFacing = TRUE,
                  adj = c(0, 0.5),
                  cex = 0.25,
                  col = ifelse(tip_tax[i] == "Devosia", "#2E86AB", "#A23B72"))
    }
  }
)

# ============================================================
# Tracks 2-21: One per BGC class (presence/absence)
# ============================================================
cat("Drawing 20 class rings...\n")

for (i in seq_along(class_order)) {
  cls <- class_order[i]
  presence <- as.integer(class_wide[[cls]] > 0)
  color_here <- class_colors[cls]

  circos.track(
    ylim = c(0, 1),
    track.height = 0.020,
    bg.border = "grey95",
    bg.col = "grey95",
    panel.fun = function(x, y) {
      # For each genome position, if present, draw a colored rectangle
      for (j in seq_len(n_tips)) {
        if (presence[j] == 1) {
          circos.rect(j - 1, 0, j, 1,
                      col = color_here,
                      border = NA)
        }
      }
      # Add class label near the start of the ring
      circos.text(0, 0.5, cls,
                  facing = "bending.inside",
                  niceFacing = TRUE,
                  adj = c(1.2, 0.5),
                  cex = 0.35,
                  col = "black")
    }
  )
}

# ============================================================
# Legend
# ============================================================
cat("Drawing legend...\n")

legend("right",
       legend = class_order,
       fill = class_colors,
       border = NA,
       cex = 0.6,
       bty = "n",
       title = "BGC class",
       title.adj = 0,
       ncol = 2)

# ============================================================
# Close
# ============================================================
circos.clear()
dev.off()
dev.off()

cat("=== [20] Done ===\n")
