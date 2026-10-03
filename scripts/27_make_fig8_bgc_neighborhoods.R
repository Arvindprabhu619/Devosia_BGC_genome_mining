#!/usr/bin/env Rscript
# Fig 8 — Bulletproof version (avoids S7 conflict)

# Clear everything
rm(list = ls())

# Load ggplot2 FIRST, before any other package
library(ggplot2)

# Now load other packages
library(dplyr)

base_dir <- "/dgxb_home/se26plsc001/arvind/Devosia_BGC"
out_dir  <- file.path(base_dir, "figures/Figure8")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gbk_files <- list(
  "C1" = file.path(base_dir,
    "results/06_antismash/raw_output/GCF_000735585.1/NZ_JNNO01000044.1.region001.gbk"),
  "C2" = file.path(base_dir,
    "results/06_antismash/raw_output/GCF_003056345.1/NZ_PYWB01000014.1.region001.gbk"),
  "C5" = file.path(base_dir,
    "results/06_antismash/raw_output/GCF_021729885.1/NZ_CP063401.1.region003.gbk")
)

titles <- c(
  C1 = "C1  RiPP-like + hydrogen-cyanide (GCF_000735585.1)",
  C2 = "C2  NRPS-like + T1PKS (GCF_003056345.1)",
  C5 = "C5  lassopeptide + RRE-containing (GCF_021729885.1)"
)

# Simple categorize using ifelse (no dplyr::case_when to avoid S7)
categorize <- function(x) {
  x <- tolower(x)
  out <- rep("other", length(x))
  out[grepl("hypothetical|unknown|duf", x) | x == ""] <- "hypothetical"
  out[grepl("nrps|nrp|synthetase|adenylation|condensation|amp-binding", x)] <- "NRPS"
  out[grepl("pks|ketosynthase|ketoreductase|acyltransferase|^ks|^at", x)] <- "PKS"
  out[grepl("lanthipeptide|lasso|ripp|precursor|thioamitide|rre", x)] <- "RiPP"
  out[grepl("cyano|hcn|hydrogen|cyanide", x)] <- "cyanide"
  out[grepl("transporter|efflux|abc|mfs|permease", x)] <- "transport"
  out[grepl("regulator|response|luxr|tetr|sigma", x)] <- "regulation"
  out
}

cat_colors <- c(
  NRPS = "#FF7F0E", PKS = "#2CA02C", RiPP = "#1F77B4",
  cyanide = "#AEC7E8", transport = "#9467BD", regulation = "#8C564B",
  hypothetical = "#D3D3D3", other = "#7F7F7F"
)

# Parse all 3 into one long data frame
all_records <- list()

for (lbl in names(gbk_files)) {
  cat("\n=== Parsing", lbl, "===\n")
  path <- gbk_files[[lbl]]
  if (!file.exists(path)) { cat("  Missing\n"); next }

  lines <- readLines(path, warn = FALSE)
  cds_idx <- grep("^\\s{5}CDS\\s", lines)
  cat("  CDS:", length(cds_idx), "\n")
  if (length(cds_idx) == 0) next

  for (i in cds_idx) {
    header <- lines[i]
    strand <- if (grepl("complement", header)) "-" else "+"
    coords <- regmatches(header, regexpr("[0-9]+\\.\\.[0-9]+", header))
    if (length(coords) == 0) next
    sp <- strsplit(coords, "\\.\\.")[[1]]
    s <- suppressWarnings(as.numeric(sp[1]))
    e <- suppressWarnings(as.numeric(sp[2]))
    if (is.na(s) || is.na(e)) next

    end_idx <- min(i + 15, length(lines))
    block <- lines[(i + 1):end_idx]
    nf <- grep("^\\s{5}\\S", block)
    if (length(nf) > 0) block <- block[1:(nf[1] - 1)]

    gene <- NA_character_; prod <- NA_character_
    for (b in block) {
      if (grepl("/gene=", b))    gene <- gsub('.*/gene="([^"]+)".*', "\\1", b)
      if (grepl("/product=", b)) prod <- gsub('.*/product="([^"]+)".*', "\\1", b)
    }
    gname <- if (!is.na(gene)) gene else if (!is.na(prod)) prod else "hypothetical"

    all_records[[length(all_records) + 1]] <- data.frame(
      candidate  = lbl,
      start      = s,
      end        = e,
      gene_label = substr(gname, 1, 25),
      stringsAsFactors = FALSE
    )
  }
}

df <- do.call(rbind, all_records)
df$start <- as.numeric(df$start)
df$end   <- as.numeric(df$end)
df$category <- categorize(df$gene_label)

cat("\nTotal genes parsed:", nrow(df), "\n")
print(table(df$candidate))
cat("start class:", class(df$start), "  end class:", class(df$end), "\n")

# --- Save 3 separate plots ---
for (lbl in names(gbk_files)) {
  sub <- df[df$candidate == lbl, , drop = FALSE]
  if (nrow(sub) == 0) next

  x_min <- min(sub$start, na.rm = TRUE)
  x_max <- max(sub$end,   na.rm = TRUE)

  p <- ggplot() +
    geom_rect(data = sub,
              aes(xmin = start, xmax = end,
                  ymin = 0.4, ymax = 0.6,
                  fill = category),
              color = "black", linewidth = 0.3) +
    scale_fill_manual(values = cat_colors) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      labels = function(x) sprintf("%.0f kb", x / 1000)
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = titles[lbl], x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid      = element_blank(),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      plot.title      = element_text(face = "bold", size = 11, hjust = 0),
      legend.position = "none"
    )

  out_png <- file.path(out_dir, paste0("Figure8_", lbl, ".png"))
  out_pdf <- file.path(out_dir, paste0("Figure8_", lbl, ".pdf"))
  ggsave(out_png, p, width = 10, height = 2, dpi = 300)
  ggsave(out_pdf, p, width = 10, height = 2)
  cat("  Saved:", out_png, "\n")
}

cat("\n=== DONE ===\n")
cat("Check:", out_dir, "\n")
