#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/21_make_fig2_definitive.R
# ============================================================
# Fig 2 — Streptomyces-style circular phylogeny. Pure circlize.
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
tree$edge.length[tree$edge.length == 0] <- 1e-7
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

# Tip taxonomy colors
tax_map <- setNames(tax$status, tax$accession)
tip_tax <- tax_map[tip_order]
tip_tax[is.na(tip_tax)] <- "Unknown"
tip_colors <- ifelse(tip_tax == "Devosia", "#2E86AB",
              ifelse(tip_tax == "Devosia_A", "#A23B72", "grey70"))

# ---------- Open devices ----------
pdf(file.path(fig_dir, "Fig2_definitive.pdf"), width = 16, height = 16)
png(file.path(fig_dir, "Fig2_definitive.png"), width = 16, height = 16,
    units = "in", res = 300)

par(mar = c(1, 1, 1, 1))

circos.par(
  start.degree = 90,
  gap.degree = 0.5,
  track.margin = c(0.003, 0.003),
  points.overflow.warning = FALSE,
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = "genome", xlim = c(0, n_tips))

# ---------- Track 1: Tip labels ----------
cat("Track 1: Labels...\n")
circos.track(
  ylim = c(0, 1),
  track.height = 0.025,
  bg.border = NA,
  panel.fun = function(x, y) {
    circos.text(
      x = seq(0.5, n_tips - 0.5),
      y = rep(0.5, n_tips),
      labels = tip_order,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.28,
      col = tip_colors
    )
  }
)

# ---------- Tracks 2..N: One ring per class ----------
cat("Adding class rings...\n")
for (i in seq_along(class_order)) {
  cls <- class_order[i]
  presence <- as.integer(class_wide[[cls]] > 0)
  color_here <- class_colors[cls]

  circos.track(
    ylim = c(0, 1),
    track.height = 0.024,
    bg.col = "#F5F5F5",
    bg.border = "white",
    panel.fun = function(x, y) {
      for (j in seq_len(n_tips)) {
        if (presence[j] == 1) {
          circos.rect(j - 1, 0, j, 1, col = color_here, border = NA)
        }
      }
      # Class label at start of ring
      circos.text(-0.5, 0.5, sprintf("%d", i),
                  facing = "bending.inside",
                  cex = 0.4, adj = c(1.1, 0.5), col = "black")
    }
  )
}

# ---------- Legend ----------
cat("Legend...\n")
circos.clear()

# Add legend on right margin
par(fig = c(0.75, 1.0, 0.1, 0.9), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()

# Class legend
legend("topright",
       legend = sprintf("%d. %s", seq_along(class_order), class_order),
       fill = class_colors,
       cex = 0.55,
       bty = "n",
       title = "BGC class",
       y.intersp = 1.15)

# Genus legend
legend("bottomright",
       legend = c("Devosia (n=101)", "Devosia_A (n=23)"),
       fill = c("#2E86AB", "#A23B72"),
       cex = 0.7,
       bty = "n",
       title = "Genus")

dev.off()
dev.off()

cat("=== [21] Done ===\n")
cat("Saved: figures/Fig2_definitive.pdf + .png\n")
