#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/21_make_fig2_definitive.R
# ============================================================
# Fig 2 — Streptomyces-style circular phylogeny with 20 class rings.
# Key: draw rectangles for ALL genomes (present=color, absent=gray).
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(circlize)
  library(dplyr)
  library(tidyr)
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
tip_colors <- ifelse(tip_tax == "Devosia", "#2E86AB",
              ifelse(tip_tax == "Devosia_A", "#A23B72", "grey70"))

# ============================================================
# Open devices
# ============================================================
pdf(file.path(fig_dir, "Fig2_definitive.pdf"), width = 14, height = 14)
png(file.path(fig_dir, "Fig2_definitive.png"), width = 14, height = 14,
    units = "in", res = 300)

par(mar = c(1, 1, 1, 1), xpd = NA)

circos.par(
  start.degree = 90,
  gap.degree = 0.3,
  track.margin = c(0.001, 0.001),
  points.overflow.warning = FALSE,
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = "g", xlim = c(0, n_tips))

# ============================================================
# Track 1: INNERMOST — Tree tips as colored bars (genus)
# This gives a visible base ring
# ============================================================
cat("Track 1: Genus strip...\n")
circos.track(
  ylim = c(0, 1),
  track.height = 0.04,
  bg.border = NA,
  panel.fun = function(x, y) {
    for (i in seq_len(n_tips)) {
      circos.rect(i - 1, 0, i, 1,
                  col = tip_colors[i], border = NA)
    }
  }
)

# ============================================================
# Tracks 2..N+1: One ring per BGC class (ALL genomes drawn)
# Present = class color, Absent = light gray
# ============================================================
cat("Adding", n_classes, "class rings (ALL genomes)...\n")
ABSENT_COLOR <- "#F0F0F0"

for (i in seq_along(class_order)) {
  cls <- class_order[i]
  presence <- as.integer(class_wide[[cls]] > 0)
  color_present <- class_colors[cls]
  n_present <- sum(presence)
  cat(sprintf("  Ring %d (%s): %d/%d\n", i, cls, n_present, n_tips))

  # Pass color_present into panel.fun environment
  local_color_present <- color_present
  circos.track(
    ylim = c(0, 1),
    track.height = 0.024,
    bg.border = NA,
    panel.fun = function(x, y) {
      for (j in seq_len(n_tips)) {
        col_j <- if (presence[j] == 1) local_color_present else ABSENT_COLOR
        circos.rect(j - 1, 0, j, 1, col = col_j, border = NA)
      }
    }
  )
}

# ============================================================
# Track (outermost): Class names as small text
# ============================================================
cat("Adding class name labels...\n")
circos.track(
  ylim = c(0, n_classes),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y) {
    for (i in seq_along(class_order)) {
      circos.text(0, i - 0.5, sprintf("%d. %s", i, class_order[i]),
                  facing = "bending.inside",
                  cex = 0.35,
                  adj = c(0, 0.5),
                  col = "black")
    }
  }
)

circos.clear()

# ============================================================
# Legend on right side (clean, not overlapping)
# ============================================================
par(fig = c(0.72, 1.0, 0.05, 0.95), new = TRUE, mar = c(2, 1, 2, 1))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, n_classes + 3))

# Title
text(0, n_classes + 2, "BGC class", adj = 0, font = 2, cex = 1.2)

# Class entries
for (i in seq_len(n_classes)) {
  y_pos <- n_classes - i + 1
  rect(0.0, y_pos - 0.35, 0.08, y_pos + 0.35,
       col = class_colors[class_order[i]], border = "grey30")
  text(0.12, y_pos, sprintf("%d. %s", i, class_order[i]),
       adj = 0, cex = 0.7)
}

# Genus section
text(0, -0.5, "Genus", adj = 0, font = 2, cex = 1.1)
rect(0.0, -1.4, 0.08, -0.8, col = "#2E86AB", border = "grey30")
text(0.12, -1.1, "Devosia (n=101)", adj = 0, cex = 0.75)
rect(0.0, -2.3, 0.08, -1.7, col = "#A23B72", border = "grey30")
text(0.12, -2.0, "Devosia_A (n=23)", adj = 0, cex = 0.75)

dev.off()
dev.off()

cat("=== [21] Done ===\n")
cat("Saved: figures/Fig2_definitive.pdf + .png\n")
