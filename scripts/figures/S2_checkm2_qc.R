#!/usr/bin/env Rscript
# =============================================================================
# S2_checkm2_qc.R — Supplementary: CheckM2 QC + completeness distribution
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
root     <- base_dir
fig_dir  <- file.path(base_dir, "figures")

PAL <- list(primary="#2166AC", secondary="#B2182B", accent="#EF8A62")
theme_pub <- function(base_size=11) {
  theme_minimal(base_size=base_size) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=base_size+1))
}

checkm <- read.delim(file.path(root, "results/03_checkm2/quality_report.tsv"),
                     stringsAsFactors = FALSE)

y_top  <- max(checkm$contamination, na.rm = TRUE)
x_left <- min(checkm$completeness, na.rm = TRUE)

pS2a <- ggplot(checkm, aes(x = completeness, y = contamination)) +
  geom_point(alpha = 0.6, color = PAL$primary, size = 2) +
  geom_hline(yintercept = 5, linetype = "dashed", color = PAL$accent) +
  geom_vline(xintercept = 95, linetype = "dashed", color = PAL$accent) +
  annotate("text", x = 95, y = y_top, label = "Completeness \u2265 95%",
           hjust = 1.1, vjust = 1.2, color = PAL$accent, size = 3.2, fontface = "bold") +
  annotate("text", x = x_left, y = 5, label = "Contamination \u2264 5%",
           hjust = -0.1, vjust = -0.6, color = PAL$accent, size = 3.2, fontface = "bold") +
  labs(x = "Completeness (%)", y = "Contamination (%)",
       title = "CheckM2 quality assessment") +
  theme_pub()

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
cat("S2 saved\n")
