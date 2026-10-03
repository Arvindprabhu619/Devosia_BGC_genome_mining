#!/usr/bin/env Rscript
# =============================================================================
# 26_make_bgc_lengths.R
# Supplementary: BGC region length distribution
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
in_file  <- file.path(base_dir, "results/06_antismash/parsed/bgc_regions.tsv")
out_dir  <- file.path(base_dir, "figures/Supplementary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

regions <- read.delim(in_file, stringsAsFactors = FALSE)
regions$length_kb <- regions$length / 1000

cat("BGCs:", nrow(regions), "\n")
cat("Length summary (kb):\n")
print(summary(regions$length_kb))

med <- median(regions$length_kb, na.rm = TRUE)

p <- ggplot(regions, aes(x = length_kb)) +
  geom_histogram(bins = 40, fill = "#67A9CF", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = med, color = "#B2182B",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = med, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Median = %.1f kb", med),
           color = "#B2182B", size = 3.5) +
  labs(
    title = "BGC region length distribution",
    x = "BGC length (kb)",
    y = "Number of BGCs",
    caption = sprintf("n = %d regions", nrow(regions))
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "S_bgc_lengths.pdf"), p, width = 8, height = 5)
ggsave(file.path(out_dir, "S_bgc_lengths.png"), p, width = 8, height = 5, dpi = 300)
cat("Saved:", out_dir, "\n")
