#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/21_make_fig2_definitive.R
# ============================================================
# Fig 2 — Circular phylogeny with 20 BGC class rings.
# Uses circlize with proper per-genome rectangle drawing.
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

# ---------- Class matrix (binary presence) ----------
# Each row = genome, each column = class
presence_matrix <- matrix(0, nrow = n_tips, ncol = 0)
class_names <- sort(unique(classes$product_class))

# Build matrix efficiently
class_wide <- classes %>%
  count(genome_accession, product_class, name = "n") %>%
  pivot_wider(id_cols = genome_accession, names_from = product_class,
              values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(class_wide) <- class_wide$genome_accession
class_wide$genome_accession <- NULL
class_wide <- class_wide[tip_order, , drop = FALSE]

# Order classes by total presence
class_order <- names(sort(colSums(class_wide > 0), decreasing = TRUE))
class_wide <- class_wide[, class_order, drop = FALSE]
n_classes <- length(class_order)
cat(sprintf("Classes: %d\n", n_classes))

# Verify data
cat(sprintf("Terpene presence: %d/%d genomes\n",
            sum(class_wide[["terpene"]] > 0), n_tips))
cat(sprintf("Hserlactone presence: %d/%d genomes\n",
            sum(class_wide[["hserlactone"]] > 0), n_tips))

# ---------- Colors ----------
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
tip_colors <- ifelse(tip_tax == "Devosia", "#2E86AB",
              ifelse(tip_tax == "Devosia_A", "#A23B72", "grey70"))

# ---------- Open devices ----------
pdf(file.path(fig_dir, "Fig2_definitive.pdf"), width = 16, height = 16)
png(file.path(fig_dir, "Fig2_definitive.png"), width = 16, height = 16,
    units = "in", res = 300)

par(mar = c(2, 2, 3, 2), xpd = NA)

circos.par(
  start.degree = 90,
  gap.degree = 0.5,
  track.margin = c(0.005, 0.005),
  points.overflow.warning = FALSE,
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = "genome", xlim = c(0, n_tips))

# ============================================================
# Track 1: Inner ring — Tip labels
# ============================================================
cat("Track 1: Tip labels...\n")
circos.track(
  ylim = c(0, 1),
  track.height = 0.04,
  bg.border = NA,
  panel.fun = function(x, y) {
    # Place tip labels around circle
    for (i in seq_len(n_tips)) {
      circos.text(
        x = i - 0.5,
        y = 0.5,
        labels = tip_order[i],
        facing = "bending.inside",
        niceFacing = TRUE,
        cex = 0.30,
        col = tip_colors[i]
      )
    }
  }
)

# ============================================================
# Tracks 2..N+1: One ring per BGC class
# ============================================================
cat("Adding", n_classes, "class rings...\n")
for (i in seq_along(class_order)) {
  cls <- class_order[i]
  presence <- as.integer(class_wide[[cls]] > 0)
  color_here <- class_colors[cls]
  n_present <- sum(presence)
  cat(sprintf("  Ring %d (%s): %d genomes present\n", i, cls, n_present))

  circos.track(
    ylim = c(0, 1),
    track.height = 0.025,
    bg.col = "#EEEEEE",
    bg.border = "white",
    panel.fun = function(x, y) {
      # In panel.fun, x is the sector's xlim range
      # We can draw rectangles using absolute positions
      for (j in seq_len(n_tips)) {
        if (presence[j] == 1) {
          circos.rect(
            xleft = j - 1,
            ybottom = 0,
            xright = j,
            ytop = 1,
            col = color_here,
            border = NA
          )
        }
      }
    }
  )
}

circos.clear()

# ============================================================
# Legend outside circos
# ============================================================
cat("Legend...\n")

# Save plot region for legend
par(fig = c(0.72, 1.0, 0.1, 0.9), new = TRUE, mar = c(0, 0, 2, 2))

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, n_classes + 2))

# Title
text(0, n_classes + 1.5, "BGC class", adj = 0, font = 2, cex = 1.1)

# Class legend
for (i in seq_len(n_classes)) {
  y_pos <- n_classes - i + 1
  rect(0, y_pos - 0.4, 0.15, y_pos + 0.4,
       col = class_colors[class_order[i]], border = "grey30")
  text(0.2, y_pos, sprintf("%d. %s", i, class_order[i]),
       adj = 0, cex = 0.75)
}

# Genus legend below
text(0, -0.3, "Genus", adj = 0, font = 2, cex = 1.0)
rect(0, -1.2, 0.15, -0.8, col = "#2E86AB", border = "grey30")
text(0.2, -1.0, "Devosia (n=101)", adj = 0, cex = 0.75)
rect(0, -2.0, 0.15, -1.6, col = "#A23B72", border = "grey30")
text(0.2, -1.8, "Devosia_A (n=23)", adj = 0, cex = 0.75)

# Close
dev.off()
dev.off()

cat("=== [21] Done ===\n")
cat("Saved: figures/Fig2_definitive.pdf + .png\n")
