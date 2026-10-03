#!/usr/bin/env Rscript
# =============================================================================
# 04_gcf_diversity.R — Figure 4: GCF diversity (3 panels)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
root     <- base_dir
fig_dir  <- file.path(base_dir, "figures")

PAL <- list(primary="#2166AC", secondary="#B2182B", accent="#EF8A62",
            core="#2166AC", rare="#B2182B", grid="grey90")
theme_pub <- function(base_size=11) {
  theme_minimal(base_size=base_size) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=base_size+1))
}

gcf_files <- list.files(file.path(root, "results/09_bigscape/parsed"),
                        pattern = "^gcf_prevalence_c.*\\.tsv$", full.names = TRUE)
gcf_list <- lapply(gcf_files, function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$cutoff <- gsub(".*c(\\d\\.\\d)\\.tsv$", "\\1", f)
  d
})
gcf_all <- bind_rows(gcf_list)
gcf_all$cutoff <- factor(gcf_all$cutoff, levels = c("0.3", "0.5", "0.7"))

# Panel A
p4a <- ggplot(gcf_all, aes(x = prevalence, fill = cutoff)) +
  geom_histogram(bins = 30, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 10, linetype = "dashed", color = PAL$rare, linewidth = 0.5) +
  geom_vline(xintercept = 90, linetype = "dashed", color = PAL$core, linewidth = 0.5) +
  facet_wrap(~ cutoff, ncol = 3,
             labeller = labeller(cutoff = function(x) paste0("c", x))) +
  scale_fill_manual(values = c("0.3"="#E8B84C", "0.5"="#B15F4A", "0.7"=PAL$primary),
                    guide = "none") +
  labs(x = "GCF prevalence (% of 124 genomes)", y = "Number of GCFs",
       title = "GCF prevalence distribution across cutoffs") +
  theme_pub()

# Panel B
gcf07 <- gcf_all %>% filter(cutoff == "0.7") %>% arrange(desc(prevalence)) %>%
  mutate(rank = row_number())
p4b <- ggplot(gcf07, aes(x = rank, y = prevalence, fill = category)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("core"=PAL$core, "intermediate"=PAL$secondary,
                                "rare"=PAL$rare), name = "Category") +
  labs(x = "GCFs (ranked by prevalence)", y = "Prevalence (%)",
       title = "GCFs ranked by prevalence (c0.7)") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel C
dev_venn <- read.delim(file.path(root, "results/09_bigscape/parsed/gcf_overlap_devosia_A.tsv"),
                       stringsAsFactors = FALSE)
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
  scale_fill_manual(values = c("Devosia\nspecific"=PAL$primary,
                                "Devosia_A\nspecific"=PAL$secondary,
                                "Shared"=PAL$accent)) +
  labs(x = NULL, y = "Number of GCFs", title = "Devosia vs Devosia_A (c0.7)") +
  theme_pub() +
  expand_limits(y = max(venn_data$n) * 1.15)

fig4 <- (p4a | p4b) / p4c +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "Fig4_gcf_diversity.pdf"), fig4,
       width = 12, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "Fig4_gcf_diversity.png"), fig4,
       width = 12, height = 8, dpi = 300)
cat("Fig 4 saved\n")
