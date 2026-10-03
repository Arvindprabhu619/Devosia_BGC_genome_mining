#!/dgxb_home/se26plsc001/miniconda3/envs/gtdbtk/bin/Rscript
# ============================================================
# scripts/13_make_figures.R
# ============================================================
# Purpose: Publication-grade manuscript figures.
#
# Outputs:
#   Fig1_workflow.pdf/png
#   Fig2_phylogenetic_landscape.pdf/png
#   Fig3_biosynthetic_diversity.pdf/png
#   Fig4_gcf_diversity.pdf/png
#   Fig5_phylogenetic_structure.pdf/png
#   Fig6_ecological_evolutionary.pdf/png
#   S1_multitrack_tree.pdf/png
#   S2_assembly_stats.pdf/png
#   S3_candidate_details.pdf/png
#   S4_devosia_A_sensitivity.pdf/png
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(cowplot); library(viridis)
  library(RColorBrewer); library(scales); library(ape)
  library(ggtree); library(ggtreeExtra); library(ggnewscale)
  library(ggpubr)
})

# ---------- Palette ----------
PAL <- list(
  primary   = "#2E86AB",
  secondary = "#A23B72",
  accent    = "#F18F01",
  rare      = "#C73E1D",
  core      = "#3B8EA5",
  neutral   = "#333333",
  grid      = "#EEEEEE",
  bg        = "#FFFFFF"
)

# ---------- Global theme ----------
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle   = element_text(size = base_size, color = "#555555"),
      plot.tag        = element_text(face = "bold", size = base_size + 4),
      axis.title      = element_text(size = base_size, face = "bold"),
      axis.text       = element_text(size = base_size - 1, color = PAL$neutral),
      legend.title    = element_text(size = base_size - 1, face = "bold"),
      legend.text     = element_text(size = base_size - 2),
      legend.key.size = unit(0.35, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = PAL$grid, size = 0.3),
      strip.text      = element_text(face = "bold", size = base_size),
      plot.margin     = margin(6, 6, 6, 6),
      plot.background = element_rect(fill = PAL$bg, color = NA)
    )
}

# ---------- Paths ----------
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(root, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== [13] Publication-grade figures ===\n")

# ---------- Load data ----------
regions  <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_regions.tsv"), stringsAsFactors = FALSE)
classes  <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_class_assignments.tsv"), stringsAsFactors = FALSE)
novelty  <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv"), stringsAsFactors = FALSE)
genomes  <- readLines(file.path(root, "results/05_gtdbtk/final_genome_list.txt"))
genomes  <- genomes[genomes != ""]
tax      <- read.delim(file.path(root, "results/05_gtdbtk/taxonomy_assignments.tsv"), stringsAsFactors = FALSE)
tree     <- read.tree(file.path(root, "results/05b_phylogeny/devosia_pruned.tree"))

# ============================================================
# FIGURE 2 — Phylogenetic landscape (BGC class ring heatmap)
# ============================================================
cat("Fig 2: Phylogenetic landscape...\n")

# Class-by-genome presence matrix
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

# Order classes by total presence
class_order <- names(sort(colSums(class_wide > 0), decreasing = TRUE))
class_wide <- class_wide[, class_order, drop = FALSE]

# Class palette
class_pal <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(ncol(class_wide)),
  class_order
)

# Per-genome BGC abundance
bgc_counts <- regions %>% count(genome_accession, name = "n_bgcs")
bgc_counts <- data.frame(
  genome_accession = tree$tip.label,
  n_bgcs = bgc_counts$n_bgcs[match(tree$tip.label, bgc_counts$genome_accession)]
)
bgc_counts$n_bgcs[is.na(bgc_counts$n_bgcs)] <- 0

# Per-genome novelty proportion
gnovel <- novelty %>%
  group_by(genome_accession) %>%
  summarise(prop_novel = sum(novelty_category == "putatively_novel") / n(),
            .groups = "drop")
gnovel <- data.frame(
  genome_accession = tree$tip.label,
  prop_novel = gnovel$prop_novel[match(tree$tip.label, gnovel$genome_accession)]
)
gnovel$prop_novel[is.na(gnovel$prop_novel)] <- 0

# Taxonomy
tax_map <- setNames(tax$status, tax$accession)
tax_status <- tax_map[tree$tip.label]
tax_status[is.na(tax_status)] <- "Unknown"
tax_df <- data.frame(genome_accession = tree$tip.label, genus = tax_status)

# Build ggtree base
p2_tree <- ggtree(tree, layout = "circular", size = 0.4, open.angle = 10) %<+%
  tax_df

# --- Build 21 rings (one per class) ---
class_long <- class_wide %>%
  tibble::rownames_to_column("genome_accession") %>%
  pivot_longer(-genome_accession, names_to = "class", values_to = "count") %>%
  mutate(present = as.integer(count > 0))

# Use discrete palette
ring_palette <- class_pal

p2 <- p2_tree
for (i in seq_along(class_order)) {
  cls <- class_order[i]
  sub <- class_long %>% filter(class == cls)
  offset_val <- 0.02 + (i - 1) * 0.025
  p2 <- p2 + geom_fruit(
    data = sub,
    geom = geom_tile,
    mapping = aes(y = genome_accession, x = "present", fill = factor(present)),
    offset = offset_val,
    pwidth = 0.022,
    color = "white", size = 0.05,
    show.legend = FALSE
  )
}

# Apply fill scale (only present shown in black; absent white)
p2 <- p2 +
  scale_fill_manual(values = c("0" = "white", "1" = "black"), guide = "none")

ggsave(file.path(fig_dir, "Fig2_phylogenetic_landscape.pdf"), p2,
       width = 14, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_phylogenetic_landscape.png"), p2,
       width = 14, height = 14, dpi = 300)

cat("  Fig 2 saved\n")

# ============================================================
# FIGURE 1 — Workflow (simplified schematic, no plot)
# ============================================================
cat("Fig 1: Workflow...\n")

workflow_data <- data.frame(
  x = rep(1, 6),
  y = 6:1,
  stage = c("516 NCBI Devosia records", "174 after curation",
            "138 passed CheckM2", "124 taxonomically validated",
            "663 BGC regions predicted", "7 priority candidates"),
  n = c(516, 174, 138, 124, 663, 7)
)

p1 <- ggplot(workflow_data, aes(x = x, y = y)) +
  geom_label(aes(label = stage), fill = PAL$primary, color = "white",
             size = 3.5, label.padding = unit(0.4, "lines"),
             label.r = unit(0.3, "lines"), fontface = "bold") +
  geom_segment(aes(x = x, xend = x, y = y - 0.35, yend = y - 0.65),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = PAL$accent, size = 1) +
  xlim(0.5, 1.5) + ylim(0.5, 6.5) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(fig_dir, "Fig1_workflow.pdf"), p1,
       width = 6, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "Fig1_workflow.png"), p1,
       width = 6, height = 8, dpi = 300)

cat("  Fig 1 saved\n")

cat("=== [13] Part 1 complete ===\n")

# ============================================================
# FIGURE 3 — Biosynthetic diversity (3 panels)
# ============================================================
cat("Fig 3: Biosynthetic diversity...\n")

# --- Panel A: Class composition ---
class_counts <- classes %>% count(product_class, sort = TRUE) %>%
  mutate(product_class = factor(product_class, levels = rev(product_class)))

p3a <- ggplot(class_counts, aes(x = n, y = product_class)) +
  geom_col(fill = PAL$primary, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.3, size = 3) +
  labs(x = "Number of class assignments", y = NULL,
       title = "Biosynthetic class composition") +
  theme_pub() +
  xlim(0, max(class_counts$n) * 1.15)

# --- Panel B: Genome × class heatmap ---
genome_class <- classes %>%
  count(genome_accession, product_class) %>%
  complete(genome_accession = genomes,
           product_class = unique(classes$product_class),
           fill = list(n = 0))

# Order classes and genomes
gm_order <- genome_class %>%
  group_by(genome_accession) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(genome_accession)

cls_order <- class_counts$product_class
genome_class$genome_accession <- factor(genome_class$genome_accession,
                                         levels = gm_order)
genome_class$product_class <- factor(genome_class$product_class,
                                      levels = cls_order)

p3b <- ggplot(genome_class, aes(x = product_class, y = genome_accession, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Count",
                       breaks = c(1, 3, 5)) +
  labs(x = NULL, y = NULL,
       title = "Genome × class heatmap (124 genomes × 21 classes)") +
  theme_pub(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

# --- Panel C: Novelty classification ---
novelty_sum <- novelty %>% count(novelty_category) %>%
  mutate(pct = 100 * n / sum(n),
         category_label = c("Known", "Putatively novel", "Related")[
           match(novelty_category, c("known", "putatively_novel", "related"))])

p3c <- ggplot(novelty_sum, aes(x = reorder(novelty_category, -n), y = n,
                                fill = novelty_category)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(
    "putatively_novel" = PAL$accent,
    "related"          = PAL$secondary,
    "known"            = PAL$primary
  )) +
  scale_x_discrete(labels = c("putatively_novel" = "Putatively\nnovel",
                               "related"          = "Related",
                               "known"            = "Known")) +
  labs(x = NULL, y = "Number of BGCs",
       title = "MIBiG-based novelty classification") +
  theme_pub() +
  expand_limits(y = max(novelty_sum$n) * 1.15)

# Composite
fig3 <- (p3a / p3c | p3b) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "Fig3_biosynthetic_diversity.pdf"), fig3,
       width = 14, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "Fig3_biosynthetic_diversity.png"), fig3,
       width = 14, height = 9, dpi = 300)

cat("  Fig 3 saved\n")

# ============================================================
# FIGURE 4 — GCF diversity and prevalence
# ============================================================
cat("Fig 4: GCF diversity...\n")

# Load GCF prevalence across cutoffs
gcf_files <- list.files(file.path(root, "results/09_bigscape/parsed"),
                        pattern = "^gcf_prevalence_c.*\\.tsv$", full.names = TRUE)
gcf_list <- lapply(gcf_files, function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$cutoff <- gsub(".*c(\\d\\.\\d)\\.tsv$", "\\1", f)
  d
})
gcf_all <- bind_rows(gcf_list)
gcf_all$cutoff <- factor(gcf_all$cutoff, levels = c("0.3", "0.5", "0.7"))

# --- Panel A: Prevalence distribution ---
p4a <- ggplot(gcf_all, aes(x = prevalence, fill = cutoff)) +
  geom_histogram(bins = 30, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 10, linetype = "dashed", color = PAL$rare, linewidth = 0.5) +
  geom_vline(xintercept = 90, linetype = "dashed", color = PAL$core, linewidth = 0.5) +
  # (old vline removed)
  facet_wrap(~ cutoff, ncol = 3, labeller = labeller(cutoff = function(x) paste0("c", x))) +
  scale_fill_manual(values = c("0.3" = "#E8B84C", "0.5" = "#B15F4A", "0.7" = PAL$primary),
                    guide = "none") +
  labs(x = "GCF prevalence (% of 124 genomes)", y = "Number of GCFs",
       title = "GCF prevalence distribution across cutoffs") +
  theme_pub()

# --- Panel B: Ranked GCFs (c0.7) ---
gcf07 <- gcf_all %>% filter(cutoff == "0.7") %>% arrange(desc(prevalence))
gcf07$rank <- seq_len(nrow(gcf07))

p4b <- ggplot(gcf07, aes(x = rank, y = prevalence, fill = category)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("core" = PAL$core,
                                "intermediate" = PAL$secondary,
                                "rare" = PAL$rare),
                    name = "Category") +
  labs(x = "GCFs (ranked by prevalence)", y = "Prevalence (%)",
       title = "GCFs ranked by prevalence (c0.7)") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# --- Panel C: Devosia vs Devosia_A ---
dev_venn <- tryCatch({
  d <- read.delim(file.path(root, "results/09_bigscape/parsed/gcf_overlap_devosia_A.tsv"),
                  stringsAsFactors = FALSE)
  d
}, error = function(e) NULL)

if (!is.null(dev_venn)) {
  venn_data <- data.frame(
    category = c("Devosia\nspecific", "Devosia_A\nspecific", "Shared"),
    n = c(dev_venn$n_gcfs[dev_venn$category == "Devosia_specific"],
          dev_venn$n_gcfs[dev_venn$category == "Devosia_A_specific"],
          dev_venn$n_gcfs[dev_venn$category == "shared"])
  )
  venn_data$category <- factor(venn_data$category,
                                levels = c("Devosia\nspecific", "Shared", "Devosia_A\nspecific"))

  p4c <- ggplot(venn_data, aes(x = category, y = n, fill = category)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = n), vjust = -0.3, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("Devosia\nspecific" = PAL$primary,
                                  "Devosia_A\nspecific" = PAL$secondary,
                                  "Shared" = PAL$accent)) +
    labs(x = NULL, y = "Number of GCFs",
         title = "Devosia vs Devosia_A (c0.7)") +
    theme_pub() +
    expand_limits(y = max(venn_data$n) * 1.15)
} else {
  p4c <- ggplot() + theme_void()
}

fig4 <- (p4a | p4b) / p4c +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "Fig4_gcf_diversity.pdf"), fig4,
       width = 12, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "Fig4_gcf_diversity.png"), fig4,
       width = 12, height = 8, dpi = 300)

cat("  Fig 4 saved\n")

# ============================================================
# FIGURE 5 — Phylogenetic structure
# ============================================================
cat("Fig 5: Phylogenetic structure...\n")

lambda_file <- file.path(root, "results/11_stats/pagels_lambda.tsv")
logreg_file <- file.path(root, "results/11_stats/logistic_regression.tsv")

# --- Panel A: Pagel's λ ---
if (file.exists(lambda_file)) {
  lambda_df <- read.delim(lambda_file, stringsAsFactors = FALSE)
  lambda_df$trait_label <- c("BGC abundance", "Novelty proportion")[
    match(lambda_df$trait, c("n_bgcs", "prop_novel"))
  ]

  p5a <- ggplot(lambda_df, aes(x = trait_label, y = lambda)) +
    geom_col(fill = PAL$primary, width = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = sprintf("lambda = %.3f\nP = %s", lambda, p_value)),
              vjust = -0.3, size = 3.2) +
    ylim(0, 1.2) +
    labs(x = NULL, y = "Pagel lambda",
         title = "Phylogenetic signal of BGC traits") +
    theme_pub()
}

# --- Panel B: Logistic regression forest plot ---
if (file.exists(logreg_file)) {
  logreg <- read.delim(logreg_file, stringsAsFactors = FALSE)
  logreg$class_label <- gsub("_", " ", logreg$class)
  logreg <- logreg %>% arrange(beta) %>% mutate(class_label = factor(class_label, levels = class_label))

  # Color by significance
  logreg$sig_color <- ifelse(logreg$fdr < 0.05, PAL$secondary, "grey70")

  p5b <- ggplot(logreg, aes(x = beta, y = class_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se,
                        color = sig_color),
                    height = 0.25, size = 0.6) +
    geom_point(aes(color = sig_color), size = 2.5) +
    scale_color_identity() +
    labs(x = "Logistic regression coefficient (beta)",
         y = NULL,
         title = "Phylogenetic associations by BGC class") +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = PAL$grid, size = 0.3))

  # Add legend manually
  p5b <- p5b +
    annotate("text", x = -4.5, y = 1, label = "* FDR < 0.05",
             color = PAL$secondary, size = 3, hjust = 0) +
    annotate("text", x = -4.5, y = 1.8, label = "* FDR >= 0.05",
             color = "grey70", size = 3, hjust = 0)
}

fig5 <- p5a | p5b +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "Fig5_phylogenetic_structure.pdf"), fig5,
       width = 14, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "Fig5_phylogenetic_structure.png"), fig5,
       width = 14, height = 7, dpi = 300)

cat("  Fig 5 saved\n")

cat("=== [13] Part 2 complete ===\n")

# ============================================================
# FIGURE 6 — Ecological + evolutionary structure
# ============================================================
cat("Fig 6: Ecological + evolutionary structure...\n")

eco_file <- file.path(root, "results/11_stats/ecological_association_classes.tsv")
age_file <- file.path(root, "results/11_stats/gcf_ages.tsv")

# --- Panel A: Niche × class associations (all, sig highlighted) ---
if (file.exists(eco_file)) {
  eco <- read.delim(eco_file, stringsAsFactors = FALSE)
  eco$class_label <- gsub("_", " ", eco$class)
  eco$niche_label <- tools::toTitleCase(eco$niche)
  eco$logp <- -log10(eco$fdr)
  eco$sig <- ifelse(eco$fdr < 0.05, "FDR < 0.05", "ns")

  # Order by significance
  eco <- eco %>% arrange(desc(logp))
  class_order_eco <- unique(eco$class_label)
  niche_order_eco <- unique(eco$niche_label)

  p6a <- ggplot(eco, aes(x = niche_label, y = class_label, fill = logp)) +
    geom_tile(color = "white", size = 0.4) +
    scale_fill_viridis_c(option = "plasma", direction = 1, name = "-log10(FDR)",
                          limits = c(0, max(eco$logp) * 1.1)) +
    labs(x = NULL, y = NULL,
         title = "BGC class × niche associations") +
    theme_pub(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank())
}

# --- Panel B: GCF evolutionary ages ---
if (file.exists(age_file)) {
  age <- read.delim(age_file, stringsAsFactors = FALSE)
  age_summary <- age %>% count(age_category) %>%
    mutate(age_category = factor(age_category,
      levels = c("ancient", "old", "intermediate", "recent", "singleton")))

  p6b <- ggplot(age_summary, aes(x = age_category, y = n, fill = age_category)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c(
      "ancient"      = "#440154",
      "old"          = "#3B528B",
      "intermediate" = "#21918C",
      "recent"       = "#5EC962",
      "singleton"    = "#FDE725"
    )) +
    labs(x = "Evolutionary age category", y = "Number of GCFs",
         title = "GCF evolutionary stratification") +
    theme_pub() +
    expand_limits(y = max(age_summary$n) * 1.15)
}

fig6 <- p6a | p6b +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "Fig6_ecological_evolutionary.pdf"), fig6,
       width = 13, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "Fig6_ecological_evolutionary.png"), fig6,
       width = 13, height = 6, dpi = 300)

cat("  Fig 6 saved\n")

# ============================================================
# FIGURE S3 — Prioritized candidate BGCs (main Fig 7)
# ============================================================
cat("Fig S3 (candidates):...\n")

cand_file <- file.path(root, "results/12_prioritization/prioritized_candidates.tsv")
if (file.exists(cand_file)) {
  cand <- read.delim(cand_file, stringsAsFactors = FALSE)
  cand$class_display <- gsub(";", " + ", cand$product_classes)
  cand$short_id <- gsub("_region", " r", cand$bgc_id)
  cand <- cand %>% mutate(short_id = factor(short_id, levels = rev(short_id)))

  # Stacked bar showing score components
  if ("reasons" %in% colnames(cand)) {
    # Parse reason string into components
    score_parts <- cand %>%
      separate_longer_delim(reasons, delim = ";") %>%
      separate(reasons, into = c("component", "value"), sep = "\\+", fill = "right") %>%
      mutate(value = as.integer(value))

    p_cand <- ggplot(score_parts, aes(x = short_id, y = value, fill = component)) +
      geom_col(width = 0.7) +
      coord_flip() +
      scale_fill_brewer(palette = "Set2", name = "Score component") +
      labs(x = NULL, y = "Prioritization score",
           title = "Prioritized BGC candidates (score breakdown)") +
      theme_pub()

    ggsave(file.path(fig_dir, "S3_candidate_details.pdf"), p_cand,
           width = 10, height = 6, dpi = 300)
    ggsave(file.path(fig_dir, "S3_candidate_details.png"), p_cand,
           width = 10, height = 6, dpi = 300)
    cat("  S3 saved\n")
  }
}

# ============================================================
# FIGURE S1 — Multitrack circular tree
# ============================================================
cat("Fig S1: Multitrack tree...\n")

track_abundance <- bgc_counts
track_novelty_df <- gnovel
track_tax_df <- tax_df

pS1 <- ggtree(tree, layout = "circular", size = 0.3, open.angle = 10) %<+% track_tax_df

# Abundance bar
pS1 <- pS1 + geom_fruit(
  data = track_abundance,
  geom = geom_col,
  mapping = aes(y = genome_accession, x = n_bgcs),
  fill = PAL$primary, offset = 0.05, pwidth = 0.3
)

# Novelty heatmap
pS1 <- pS1 + geom_fruit(
  data = track_novelty_df,
  geom = geom_tile,
  mapping = aes(y = genome_accession, x = "n", fill = prop_novel),
  offset = 0.4, pwidth = 0.08
)

pS1 <- pS1 +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Novelty\nproportion")

ggsave(file.path(fig_dir, "S1_multitrack_tree.pdf"), pS1,
       width = 12, height = 12, dpi = 300)
ggsave(file.path(fig_dir, "S1_multitrack_tree.png"), pS1,
       width = 12, height = 12, dpi = 300)

cat("  S1 saved\n")

# ============================================================
# FIGURE S2 — Assembly statistics
# ============================================================
cat("Fig S2: Assembly stats...\n")

manifest <- read.delim(file.path(root, "data/metadata/genome_manifest.tsv"),
                        stringsAsFactors = FALSE)
checkm <- read.delim(file.path(root, "results/03_checkm2/quality_report.tsv"),
                      stringsAsFactors = FALSE)

if (nrow(checkm) > 0 && "completeness" %in% colnames(checkm)) {
  pS2a <- ggplot(checkm, aes(x = completeness, y = contamination)) +
    geom_point(alpha = 0.6, color = PAL$primary, size = 2) +
    geom_hline(yintercept = 5, linetype = "dashed", color = PAL$accent) +
    geom_vline(xintercept = 95, linetype = "dashed", color = PAL$accent) +
    labs(x = "Completeness (%)", y = "Contamination (%)",
         title = "CheckM2 quality assessment") +
    theme_pub()
}

if (nrow(checkm) > 0 && "completeness" %in% colnames(checkm)) {
  pS2b <- ggplot(checkm, aes(x = completeness)) +
    geom_histogram(bins = 30, fill = PAL$secondary, color = "white") +
    labs(x = "Completeness (%)", y = "Number of genomes",
         title = "Completeness distribution") +
    theme_pub()

  figS2 <- pS2a | pS2b
  ggsave(file.path(fig_dir, "S2_assembly_stats.pdf"), figS2,
         width = 11, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "S2_assembly_stats.png"), figS2,
         width = 11, height = 5, dpi = 300)
  cat("  S2 saved\n")
}

# ============================================================
# FIGURE S4 — Devosia_A sensitivity
# ============================================================
cat("Fig S4: Devosia_A sensitivity...\n")

# Novelty breakdown by taxonomy
novelty_tax <- novelty %>%
  left_join(tax %>% select(accession, status), by = c("genome_accession" = "accession")) %>%
  count(status, novelty_category) %>%
  group_by(status) %>%
  mutate(pct = 100 * n / sum(n))

pS4 <- ggplot(novelty_tax, aes(x = status, y = pct, fill = novelty_category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = c(
    "putatively_novel" = PAL$accent,
    "related"          = PAL$secondary,
    "known"            = PAL$primary
  ), name = "Novelty category") +
  labs(x = "Taxonomic assignment (GTDB R232)", y = "Proportion of BGCs (%)",
       title = "Novelty classification by genus") +
  theme_pub()

ggsave(file.path(fig_dir, "S4_devosia_A_sensitivity.pdf"), pS4,
       width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "S4_devosia_A_sensitivity.png"), pS4,
       width = 8, height = 5, dpi = 300)

cat("  S4 saved\n")

cat("\n=== [13] All figures generated ===\n")
cat("Output: ", fig_dir, "\n")
print(list.files(fig_dir, pattern = "\\.(pdf|png)$"))
