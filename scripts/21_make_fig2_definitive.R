#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/21_make_fig2_definitive.R
# ============================================================
# Fig 2 — Circular phylogeny with 20 class rings.
# FIX: thicker rings, no margin, tree in center, clean layout.
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

par(mar = c(0.5, 0.5, 0.5, 0.5), xpd = NA)

circos.par(
  start.degree = 90,
  gap.degree = 0.3,
  track.margin = c(0.0005, 0.0005),  # tiny gap
  points.overflow.warning = FALSE,
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = "g", xlim = c(0, n_tips))

# ============================================================
# Track 1: INNERMOST — Genus strip (color band)
# ============================================================
cat("Track 1: Genus strip...\n")
circos.track(
  ylim = c(0, 1),
  track.height = 0.03,
  bg.border = "white",
  panel.fun = function(x, y) {
    for (i in seq_len(n_tips)) {
      circos.rect(i - 1, 0, i, 1, col = tip_colors[i], border = NA)
    }
  }
)

# ============================================================
# Tracks 2..N+1: Class rings — ALL genomes drawn
# ============================================================
cat("Adding", n_classes, "class rings...\n")
ABSENT_COLOR <- "#E8E8E8"

# Build environment for closures
presence_list <- lapply(class_order, function(cls) as.integer(class_wide[[cls]] > 0))
names(presence_list) <- class_order

for (i in seq_along(class_order)) {
  cls <- class_order[i]
  presence <- presence_list[[cls]]
  color_present <- class_colors[cls]
  n_present <- sum(presence)
  cat(sprintf("  Ring %d (%s): %d/%d\n", i, cls, n_present, n_tips))

  circos.track(
    ylim = c(0, 1),
    track.height = 0.030,          # thicker ring
    bg.border = NA,
    panel.fun = function(x, y) {
      for (j in seq_len(n_tips)) {
        col_j <- if (presence[j] == 1) color_present else ABSENT_COLOR
        circos.rect(j - 1, 0, j, 1, col = col_j, border = NA)
      }
    }
  )
}

circos.clear()

# ============================================================
# Legend on right
# ============================================================
par(fig = c(0.75, 1.0, 0.10, 0.90), new = TRUE, mar = c(1, 1, 1, 1))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, n_classes + 3))

text(0, n_classes + 2, "BGC class", adj = 0, font = 2, cex = 1.2)

for (i in seq_len(n_classes)) {
  y_pos <- n_classes - i + 1
  rect(0.0, y_pos - 0.35, 0.08, y_pos + 0.35,
       col = class_colors[class_order[i]], border = "grey30")
  text(0.12, y_pos, sprintf("%d. %s", i, class_order[i]),
       adj = 0, cex = 0.68)
}

text(0, -0.7, "Genus", adj = 0, font = 2, cex = 1.1)
rect(0.0, -1.6, 0.08, -1.0, col = "#2E86AB", border = "grey30")
text(0.12, -1.3, "Devosia (n=101)", adj = 0, cex = 0.72)
rect(0.0, -2.5, 0.08, -1.9, col = "#A23B72", border = "grey30")
text(0.12, -2.2, "Devosia_A (n=23)", adj = 0, cex = 0.72)

# Panel title
par(fig = c(0, 1, 0.93, 1.0), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
text(0.02, 0.5, "A", adj = 0, cex = 1.8, font = 2)
text(0.06, 0.5, "Phylogenetic distribution of 20 BGC classes across 124 Devosia genomes",
     adj = 0, cex = 1.1, font = 2)

dev.off()
dev.off()

cat("=== [21] Done ===\n")
