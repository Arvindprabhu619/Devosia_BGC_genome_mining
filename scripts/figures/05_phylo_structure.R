#!/usr/bin/env Rscript
# =============================================================================
# 05_phylo_structure.R — Figure 5: Phylogenetic structure (2 panels)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
root     <- base_dir
fig_dir  <- file.path(base_dir, "figures")

PAL <- list(primary="#2166AC", secondary="#B2182B", grid="grey90")
theme_pub <- function(base_size=11) {
  theme_minimal(base_size=base_size) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=base_size+1))
}

lambda_file <- file.path(root, "results/11_stats/pagels_lambda.tsv")
logreg_file <- file.path(root, "results/11_stats/logistic_regression.tsv")

# Panel A
lambda_df <- read.delim(lambda_file, stringsAsFactors = FALSE)
lambda_df$trait_label <- c("BGC abundance", "Novelty proportion")[
  match(lambda_df$trait, c("n_bgcs", "prop_novel"))]

p5a <- ggplot(lambda_df, aes(x = trait_label, y = lambda)) +
  geom_col(fill = PAL$primary, width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = sprintf("lambda = %.3f\nP = %s", lambda, p_value)),
            vjust = -0.3, size = 3.2) +
  ylim(0, 1.2) +
  labs(x = NULL, y = "Pagel lambda",
       title = "Phylogenetic signal of BGC traits") +
  theme_pub()

# Panel B
logreg <- read.delim(logreg_file, stringsAsFactors = FALSE)
logreg$class_label <- gsub("_", " ", logreg$class)
logreg <- logreg %>% arrange(beta) %>%
  mutate(class_label = factor(class_label, levels = class_label))
logreg$sig_color <- ifelse(logreg$fdr < 0.05, PAL$secondary, "grey70")

p5b <- ggplot(logreg, aes(x = beta, y = class_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se,
                     color = sig_color),
                 height = 0.25, linewidth = 0.6) +
  geom_point(aes(color = sig_color), size = 2.5) +
  scale_color_identity() +
  labs(x = "Logistic regression coefficient (beta)", y = NULL,
       title = "Phylogenetic associations by BGC class") +
  theme_pub() +
  theme(panel.grid.major.y = element_line(color = PAL$grid, linewidth = 0.3)) +
  annotate("text", x = -4.5, y = 1, label = "* FDR < 0.05",
           color = PAL$secondary, size = 3, hjust = 0) +
  annotate("text", x = -4.5, y = 1.8, label = "* FDR >= 0.05",
           color = "grey70", size = 3, hjust = 0)

fig5 <- p5a | p5b +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "Fig5_phylogenetic_structure.pdf"), fig5,
       width = 14, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "Fig5_phylogenetic_structure.png"), fig5,
       width = 14, height = 7, dpi = 300)
cat("Fig 5 saved\n")
