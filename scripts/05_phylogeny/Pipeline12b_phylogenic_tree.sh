#!/usr/bin/env bash

# ============================================================
# Pipeline 12b — Publication-quality Devosia phylogeny
#
# Purpose:
#   Prune the five outgroup genomes from the rooted GTDB-Tk
#   tree and generate publication-quality circular phylogeny
#   outputs.
#
# Portability:
#   Repository root is inferred from this script's location.
#   No hard-coded HOME or project-specific paths.
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

PHYLO_DIR="$PROJECT_ROOT/phylogeny_110"

mkdir -p "$PHYLO_DIR"

export DEVOSIA_PROJECT_ROOT="$PROJECT_ROOT"

cat > "$PHYLO_DIR/make_publication_tree_CIRCULAR_v3.R" <<'RSCRIPT'
#!/usr/bin/env Rscript

# ============================================================
# PUBLICATION-QUALITY PHYLOGENOMIC TREE — CIRCULAR LAYOUT (v3)
# 110 Devosia genomes
#
# FIX in this version (vs v2):
# v2 depended on the 'cowplot' package, which is not installed
# in this environment (Error: there is no package called
# 'cowplot'). Rather than requiring a new package install on
# the HPC, the scale-bar inset is now composited using base
# 'grid' graphics only (grid ships with every R installation,
# no install needed). The scale bar is still drawn as an
# independent mini-plot in absolute figure-fraction coordinates
# (via grid::viewport), so it remains structurally impossible
# for it to land on any tree branch or label, regardless of
# tree shape — same guarantee as v2, without the dependency.
#
# GTDB-Tk BAC120 de novo phylogeny
# Rooted using Alphaproteobacteria outgroup
# Five outgroup genomes removed from display
# ============================================================


# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

suppressPackageStartupMessages({
    library(ape)
    library(ggtree)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(grid)   # base R package — always available, no install needed
})


# ------------------------------------------------------------
# 2. PATHS
# ------------------------------------------------------------

project_root <- Sys.getenv("DEVOSIA_PROJECT_ROOT")

if (project_root == "") {
    stop(
        "ERROR: DEVOSIA_PROJECT_ROOT environment variable is not set.\n",
        "Pipeline 12b must be executed through the repository wrapper."
    )
}

base_dir <- file.path(
    project_root,
    "phylogeny_110"
)

decorated_tree_file <- file.path(base_dir, "gtdbtk_final_rooted_v2", "gtdbtk.bac120.decorated.tree")
tax_file             <- file.path(base_dir, "Devosia_110_taxonomy.tsv")
pruned_tree_out      <- file.path(base_dir, "Devosia_110_final_rooted_pruned.tree")
metadata_file        <- file.path(base_dir, "Devosia_110_phylogeny_labels.tsv")
pdf_file             <- file.path(base_dir, "Devosia_110_phylogeny_CIRCULAR_v3.pdf")
png_file             <- file.path(base_dir, "Devosia_110_phylogeny_CIRCULAR_v3_600dpi.png")
tiff_file            <- file.path(base_dir, "Devosia_110_phylogeny_CIRCULAR_v3_600dpi.tiff")


# ------------------------------------------------------------
# 3. SETTINGS
# ------------------------------------------------------------

EXPECTED_FULL_TIPS <- 115
EXPECTED_FINAL_TIPS <- 110

SUPPORT_THRESHOLD <- 0.70

FIG_WIDTH  <- 22
FIG_HEIGHT <- 22

TIP_LABEL_SIZE <- 2.6
SUPPORT_LABEL_SIZE <- 1.7

TITLE_SIZE <- 15
SUBTITLE_SIZE <- 9.5

CHAR_WIDTH_FACTOR <- 0.017
LABEL_OFFSET_FACTOR <- 0.010

# Scale-bar inset position, in figure-fraction coordinates
# (0,0) = bottom-left corner of the whole canvas, (1,1) = top-
# right. Tune these if your longest labels ever reach into the
# bottom-left corner.
SCALE_INSET_X <- 0.03
SCALE_INSET_Y <- 0.03
SCALE_INSET_WIDTH <- 0.16
SCALE_INSET_HEIGHT <- 0.05


# ------------------------------------------------------------
# 4. OUTGROUP TIPS
# ------------------------------------------------------------

outgroup_tips <- c(
    "GCF_000743515.1_ASM74351v1_genomic",
    "GCF_001185205.1_ASM118520v1_genomic",
    "GCF_014083885.1_ASM1408388v1_genomic",
    "GCF_000689495.1_DevNy.1_genomic",
    "GCF_900119845.1_IMG-taxon_2596583700_annotated_assembly_genomic"
)


# ------------------------------------------------------------
# 5. CHECK INPUT FILES
# ------------------------------------------------------------

cat("\n============================================================\n")
cat(" Devosia 110-genome publication phylogeny (circular, v3)\n")
cat("============================================================\n\n")

if (!file.exists(decorated_tree_file)) stop("ERROR: GTDB-Tk tree not found:\n", decorated_tree_file)
if (!file.exists(tax_file)) stop("ERROR: Taxonomy file not found:\n", tax_file)

cat("Tree:\n  ", decorated_tree_file, "\n\n", sep = "")
cat("Taxonomy:\n  ", tax_file, "\n\n", sep = "")


# ------------------------------------------------------------
# 6. READ ROOTED GTDB-TK TREE
# ------------------------------------------------------------

cat("Reading GTDB-Tk tree...\n")
tree_full <- read.tree(decorated_tree_file)
cat("Input tree tips: ", length(tree_full$tip.label), "\n", sep = "")

if (length(tree_full$tip.label) != EXPECTED_FULL_TIPS) {
    stop("ERROR: Expected ", EXPECTED_FULL_TIPS, " tips but found ", length(tree_full$tip.label), ".")
}
if (!is.rooted(tree_full)) stop("ERROR: Input GTDB-Tk tree is not rooted.")
cat("Input tree: ROOTED\n\n")


# ------------------------------------------------------------
# 7. VERIFY OUTGROUPS
# ------------------------------------------------------------

cat("Checking outgroup tips...\n")
missing_outgroups <- setdiff(outgroup_tips, tree_full$tip.label)
if (length(missing_outgroups) > 0) {
    cat("\nMissing outgroup tips:\n")
    cat(paste(missing_outgroups, collapse = "\n"), "\n")
    stop("ERROR: One or more expected outgroup labels were not found.")
}
cat("All ", length(outgroup_tips), " outgroup tips found.\n\n", sep = "")


# ------------------------------------------------------------
# 8. PRUNE OUTGROUPS
# ------------------------------------------------------------

cat("Pruning outgroup tips...\n")
tree_pruned <- drop.tip(tree_full, outgroup_tips)
cat("Final tree tips: ", length(tree_pruned$tip.label), "\n", sep = "")

if (length(tree_pruned$tip.label) != EXPECTED_FINAL_TIPS) {
    stop("ERROR: Expected ", EXPECTED_FINAL_TIPS, " tips after pruning but found ", length(tree_pruned$tip.label), ".")
}
if (!is.rooted(tree_pruned)) stop("ERROR: Final tree is not recognized as rooted.")
cat("Final tree: ROOTED\n\n")


# ------------------------------------------------------------
# 9. SAVE ROOTED 110-GENOME TREE
# ------------------------------------------------------------

write.tree(tree_pruned, file = pruned_tree_out)
cat("Final rooted tree written:\n  ", pruned_tree_out, "\n\n", sep = "")


# ------------------------------------------------------------
# 10. READ TAXONOMY
# ------------------------------------------------------------

cat("Reading taxonomy metadata...\n")
taxonomy <- read.delim(tax_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("Taxonomy rows: ", nrow(taxonomy), "\n", sep = "")
cat("Taxonomy columns:\n  ", paste(colnames(taxonomy), collapse = ", "), "\n\n", sep = "")


# ------------------------------------------------------------
# 11. IDENTIFY GENOME-ID COLUMN
# ------------------------------------------------------------

possible_id_columns <- c("user_genome", "genome_id", "accession", "tip.label", "genome")
id_col <- intersect(possible_id_columns, colnames(taxonomy))
if (length(id_col) == 0) stop("ERROR: Could not identify genome-ID column.")
id_col <- id_col[1]
cat("Genome-ID column: ", id_col, "\n\n", sep = "")


# ------------------------------------------------------------
# 12. VERIFY TREE/TAXONOMY MATCH
# ------------------------------------------------------------

tree_genomes <- tree_pruned$tip.label
tax_genomes <- taxonomy[[id_col]]
missing_tax <- setdiff(tree_genomes, tax_genomes)
extra_tax <- setdiff(tax_genomes, tree_genomes)

cat("Tree genomes: ", length(tree_genomes), "\n", sep = "")
cat("Taxonomy genomes: ", length(tax_genomes), "\n", sep = "")
cat("Tree genomes missing taxonomy: ", length(missing_tax), "\n", sep = "")
cat("Taxonomy genomes absent from tree: ", length(extra_tax), "\n\n", sep = "")

if (length(missing_tax) > 0) stop("ERROR: Tree tips without taxonomy detected.")


# ------------------------------------------------------------
# 13. ALIGN TAXONOMY TO TREE ORDER
# ------------------------------------------------------------

metadata <- taxonomy[match(tree_genomes, tax_genomes), , drop = FALSE]
metadata[[id_col]] <- tree_genomes


# ------------------------------------------------------------
# 14. BUILD PUBLICATION TIP LABELS
# ------------------------------------------------------------

if ("species" %in% colnames(metadata)) {
    species_raw <- metadata$species
} else {
    species_raw <- rep(NA_character_, nrow(metadata))
}

clean_species <- gsub("^s__", "", species_raw)
clean_species[is.na(clean_species) | clean_species == "" | clean_species == "s__"] <- "Devosia sp."
clean_species <- str_trim(clean_species)

short_accession <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*$", "\\1", tree_genomes)

metadata$label_display <- paste0(short_accession, " | ", clean_species)


# ------------------------------------------------------------
# 15. CHECK LABELS
# ------------------------------------------------------------

if (anyDuplicated(metadata$label_display)) {
    duplicated_labels <- unique(metadata$label_display[duplicated(metadata$label_display)])
    stop("ERROR: Duplicate publication tip labels detected:\n", paste(duplicated_labels, collapse = "\n"))
}


# ------------------------------------------------------------
# 16. WRITE METADATA
# ------------------------------------------------------------

write.table(metadata, file = metadata_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
cat("Metadata written:\n  ", metadata_file, "\n\n", sep = "")


# ------------------------------------------------------------
# 17. TREE DEPTH / BRANCH LENGTH CHECK
# ------------------------------------------------------------

if (is.null(tree_pruned$edge.length)) stop("ERROR: Tree has no branch lengths.")
root_to_tip <- node.depth.edgelength(tree_pruned)
tree_depth <- max(root_to_tip, na.rm = TRUE)
cat("Maximum root-to-tip distance: ", signif(tree_depth, 7), " substitutions/site\n\n", sep = "")


# ------------------------------------------------------------
# 18. BUILD GGTREE — CIRCULAR LAYOUT
# ------------------------------------------------------------

cat("Building circular ggtree object...\n")

p <- ggtree(
    tree_pruned,
    layout = "circular",
    linewidth = 0.35
) %<+% metadata


# ------------------------------------------------------------
# 19. ACCESS GGTREE PLOTTING DATA
# ------------------------------------------------------------

plot_data <- p$data
required_columns <- c("node", "isTip", "x", "y", "label")
missing_columns <- setdiff(required_columns, colnames(plot_data))
if (length(missing_columns) > 0) {
    stop("ERROR: Required ggtree columns missing:\n", paste(missing_columns, collapse = ", "))
}


# ------------------------------------------------------------
# 20. EXTRACT INTERNAL NODE SUPPORT
# ------------------------------------------------------------

support_df <- plot_data %>%
    filter(isTip == FALSE) %>%
    mutate(support_numeric = suppressWarnings(as.numeric(label))) %>%
    filter(!is.na(support_numeric))

cat("Internal nodes: ", sum(plot_data$isTip == FALSE), "\n", sep = "")
cat("Numeric support values: ", nrow(support_df), "\n", sep = "")

if (nrow(support_df) == 0) {
    stop("\nERROR: No numeric internal-node support values detected.\nInspect the GTDB-Tk decorated tree.")
}


# ------------------------------------------------------------
# 21. VALIDATE SUPPORT SCALE
# ------------------------------------------------------------

support_min <- min(support_df$support_numeric, na.rm = TRUE)
support_max <- max(support_df$support_numeric, na.rm = TRUE)
cat("Support range: ", round(support_min, 4), " - ", round(support_max, 4), "\n", sep = "")

if (support_max > 1.01) {
    stop(
        "\nERROR: Support values appear to be on a 0-100 scale.\n",
        "Maximum detected value: ", support_max, "\n",
        "Do not proceed until the node-label interpretation is verified."
    )
}


# ------------------------------------------------------------
# 22. SELECT SUPPORT >= 0.70
# ------------------------------------------------------------

support_df <- support_df %>%
    filter(support_numeric >= SUPPORT_THRESHOLD) %>%
    mutate(support_display = sprintf("%.2f", support_numeric))

cat("Support values displayed (>= ", SUPPORT_THRESHOLD, "): ", nrow(support_df), "\n\n", sep = "")


# ------------------------------------------------------------
# 23. TIP LABELS
# ------------------------------------------------------------

p <- p +
    geom_tiplab2(
        aes(label = label_display),
        size = TIP_LABEL_SIZE,
        align = TRUE,
        linetype = "dotted",
        linesize = 0.18,
        offset = tree_depth * LABEL_OFFSET_FACTOR,
        hjust = 0,
        fontface = "plain"
    )


# ------------------------------------------------------------
# 24. SUPPORT LABELS
# ------------------------------------------------------------

if (nrow(support_df) > 0) {
    p <- p +
        geom_text(
            data = support_df,
            aes(x = x, y = y, label = support_display),
            size = SUPPORT_LABEL_SIZE,
            hjust = 1.15,
            vjust = -0.25,
            inherit.aes = FALSE
        )
}


# ------------------------------------------------------------
# 25. TITLE / SUBTITLE
# ------------------------------------------------------------

p <- p +
    ggtitle(
        "Phylogenomic relationships among 110 Devosia genomes",
        subtitle = "GTDB-Tk BAC120 phylogeny; rooted using an Alphaproteobacteria outgroup and displayed after outgroup pruning"
    )


# ------------------------------------------------------------
# 26. PUBLICATION THEME
# ------------------------------------------------------------

p <- p +
    theme_tree() +
    theme(
        plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0, margin = margin(b = 8)),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )


# ------------------------------------------------------------
# 27. RESERVE RADIAL LABEL SPACE
# ------------------------------------------------------------

max_label_chars <- max(nchar(metadata$label_display))
label_space <- tree_depth * CHAR_WIDTH_FACTOR * max_label_chars
x_extension <- label_space * 1.05

if (!is.finite(x_extension) || x_extension <= 0) x_extension <- tree_depth * 0.3

cat("Longest label: ", max_label_chars, " characters -> reserved radius extension: ", signif(x_extension, 4), "\n\n", sep = "")

p <- p +
    xlim(0, tree_depth + x_extension)


# ------------------------------------------------------------
# 28. BUILD SCALE-BAR MINI-PLOT (no cowplot needed)
# It is composited via grid::viewport at save time (step 29),
# in absolute figure-fraction coordinates — independent of the
# tree's polar coordinate system, so it cannot land on a branch.
# ------------------------------------------------------------

scale_width <- tree_depth * 0.15
if (!is.finite(scale_width) || scale_width <= 0) scale_width <- 0.01

scale_bar_df <- data.frame(x = c(0, scale_width), y = c(0, 0))

scale_plot <- ggplot(scale_bar_df, aes(x = x, y = y)) +
    geom_line(linewidth = 0.6) +
    annotate(
        "text",
        x = scale_width / 2,
        y = 0.45,
        label = signif(scale_width, 4),
        size = 3.2
    ) +
    annotate(
        "text",
        x = scale_width / 2,
        y = -0.45,
        label = "substitutions/site",
        size = 3.2
    ) +
    ylim(-1, 1) +
    theme_void()


# ------------------------------------------------------------
# 29. COMPOSITE + SAVE FUNCTION (grid-based, no cowplot)
# Draws the tree, then overlays the scale-bar mini-plot in a
# grid::viewport placed at absolute figure-fraction coordinates
# (SCALE_INSET_X/Y/WIDTH/HEIGHT), then closes the device.
# ------------------------------------------------------------

draw_composite <- function() {
    print(p, newpage = TRUE)
    vp <- viewport(
        x = unit(SCALE_INSET_X + SCALE_INSET_WIDTH / 2, "npc"),
        y = unit(SCALE_INSET_Y + SCALE_INSET_HEIGHT / 2, "npc"),
        width  = unit(SCALE_INSET_WIDTH, "npc"),
        height = unit(SCALE_INSET_HEIGHT, "npc")
    )
    pushViewport(vp)
    print(scale_plot, vp = vp, newpage = FALSE)
    popViewport()
}


# ------------------------------------------------------------
# 30. SAVE VECTOR PDF
# ------------------------------------------------------------

cat("Writing vector PDF...\n")
cairo_pdf(filename = pdf_file, width = FIG_WIDTH, height = FIG_HEIGHT)
draw_composite()
dev.off()


# ------------------------------------------------------------
# 31. SAVE 600 DPI PNG
# ------------------------------------------------------------

cat("Writing 600-dpi PNG...\n")
png(
    filename = png_file, width = FIG_WIDTH, height = FIG_HEIGHT,
    units = "in", res = 600, type = "cairo"
)
draw_composite()
dev.off()


# ------------------------------------------------------------
# 32. SAVE 600 DPI TIFF
# ------------------------------------------------------------

cat("Writing 600-dpi TIFF...\n")
tiff(
    filename = tiff_file, width = FIG_WIDTH, height = FIG_HEIGHT,
    units = "in", res = 600, compression = "lzw", type = "cairo"
)
draw_composite()
dev.off()


# ------------------------------------------------------------
# 33. VALIDATE OUTPUTS
# ------------------------------------------------------------

output_files <- c(pruned_tree_out, metadata_file, pdf_file, png_file, tiff_file)
missing_outputs <- output_files[!file.exists(output_files)]
if (length(missing_outputs) > 0) {
    stop("ERROR: Output files missing:\n", paste(missing_outputs, collapse = "\n"))
}


# ------------------------------------------------------------
# 34. FINAL REPORT
# ------------------------------------------------------------

cat("\n============================================================\n")
cat(" CIRCULAR PUBLICATION TREE (v3) COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Input tips:              ", length(tree_full$tip.label), "\n", sep = "")
cat("Outgroups removed:       ", length(outgroup_tips), "\n", sep = "")
cat("Final Devosia tips:      ", length(tree_pruned$tip.label), "\n", sep = "")
cat("Rooted:                  ", is.rooted(tree_pruned), "\n", sep = "")
cat("Support threshold:       >= ", SUPPORT_THRESHOLD, "\n", sep = "")
cat("Maximum distance:        ", signif(tree_depth, 7), " substitutions/site\n\n", sep = "")

cat("OUTPUT FILES\n\n")
cat("Rooted Newick:\n  ", pruned_tree_out, "\n\n", sep = "")
cat("Metadata:\n  ", metadata_file, "\n\n", sep = "")
cat("Vector PDF (use this on the slide):\n  ", pdf_file, "\n\n", sep = "")
cat("600-dpi PNG:\n  ", png_file, "\n\n", sep = "")
cat("600-dpi TIFF:\n  ", tiff_file, "\n\n", sep = "")

cat("SUCCESS: 110 / 110 genomes plotted.\n")
cat("SUCCESS: rooted tree retained.\n")
cat("SUCCESS: scale bar moved off-tree into blank corner inset (via grid, no cowplot dependency).\n")
cat("SUCCESS: circular publication-quality outputs generated.\n")
cat("\n============================================================\n")
RSCRIPT

echo
echo "============================================================"
echo "Running publication-tree R script"
echo "============================================================"
echo

Rscript "$PHYLO_DIR/make_publication_tree_CIRCULAR_v3.R"

echo
echo "============================================================"
echo "Pipeline 12b completed."
echo "============================================================"
echo
echo "Phylogeny directory:"
echo "  $PHYLO_DIR"

