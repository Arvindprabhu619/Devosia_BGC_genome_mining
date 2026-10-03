#!/usr/bin/env Rscript
# =============================================================================
# S4_novelty_by_genus.R — Supplementary: Novelty classification by genus
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
root     <- base_dir
fig_dir  <- file.path(base_dir, "figures")

PAL <- list(primary="#2166AC", secondary="#B2182B", accent="#E58601")
theme_pub <- function(base_size=11) {
  theme_minimal(base_size=base_size) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=base_size+1))
}

novelty <- read.delim(file.path(root, "results/06_antismash/parsed/bgc_novelty.tsv"),
                      stringsAsFactors = FALSE)
tax <- read.delim(file.path(root, "results/04_gtdbtk/taxonomy_assignments.tsv"),
                  stringsAsFactors = FALSE)
# tax has: accession, status
tax <- tax %>% select(accession, status)

novelty_tax <- novelty %>%
  left_join(tax, by = c("genome_accession" = "accession")) %>%
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
cat("S4 saved\n")
