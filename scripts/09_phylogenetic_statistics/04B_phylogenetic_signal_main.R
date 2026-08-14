#!/usr/bin/env Rscript

library(ape)
library(phylolm)

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS 4B: PRIMARY PHYLOGENETIC SIGNAL ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# ============================================================
# INPUTS
# ============================================================

tree_file <- "../phylogeny_110/Devosia_110_final_rooted_pruned.tree"
s2_file   <- "../manuscript_outputs/TableS2_genome_by_class_matrix.tsv"

tree <- read.tree(tree_file)

s2 <- read.delim(
    s2_file,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("\nTree tips:", length(tree$tip.label), "\n")
cat("S2 rows :", nrow(s2), "\n")

# ============================================================
# ACCESSION RECONCILIATION
# ============================================================

tree_acc <- sub(
    "^((GC[FA]_[0-9]+\\.[0-9]+)).*$",
    "\\1",
    tree$tip.label
)

stopifnot(length(tree_acc) == 110)
stopifnot(length(unique(tree_acc)) == 110)
stopifnot(length(intersect(tree_acc, s2$accession)) == 110)

s2 <- s2[match(tree_acc, s2$accession), ]

stopifnot(all(s2$accession == tree_acc))

cat("PASS: 110/110 genomes matched and aligned.\n")

# ============================================================
# ZERO-BRANCH CORRECTION
# ============================================================

EPSILON <- 1e-8

zero_before <- sum(tree$edge.length <= 0)

tree$edge.length[tree$edge.length <= 0] <- EPSILON

zero_after <- sum(tree$edge.length <= 0)

cat("\nZero-length branches before:", zero_before, "\n")
cat("Zero-length branches after :", zero_after, "\n")
cat("Minimum branch length      :", min(tree$edge.length), "\n")

# ============================================================
# BGC CLASSES
# ============================================================

classes <- setdiff(
    colnames(s2),
    "accession"
)

cat("\nTotal BGC classes:", length(classes), "\n")
print(classes)

# ============================================================
# RESULTS CONTAINER
# ============================================================

results <- list()

# ============================================================
# MAIN LOOP
# ============================================================

for (cls in classes) {

    y <- as.numeric(s2[[cls]] > 0)

    present <- sum(y == 1, na.rm = TRUE)
    absent  <- sum(y == 0, na.rm = TRUE)
    missing <- sum(is.na(y))

    cat("\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat("CLASS:", cls, "\n")
    cat("Present:", present, "\n")
    cat("Absent :", absent, "\n")
    cat("NA     :", missing, "\n")

    # --------------------------------------------------------
    # Exclude invariant / extremely rare traits
    # --------------------------------------------------------

    if (
        missing > 0 ||
        present < 5 ||
        absent < 5
    ) {

        cat("STATUS: EXCLUDED\n")

        results[[length(results) + 1]] <- data.frame(
            class = cls,
            present = present,
            absent = absent,
            n = length(y),
            estimate = NA_real_,
            SE = NA_real_,
            z = NA_real_,
            p = NA_real_,
            alpha = NA_real_,
            status = "excluded",
            stringsAsFactors = FALSE
        )

        next
    }

    # --------------------------------------------------------
    # Data frame with EXACT tree tip labels
    # --------------------------------------------------------

    dat <- data.frame(
        y = y,
        row.names = tree$tip.label
    )

    stopifnot(identical(rownames(dat), tree$tip.label))

    # --------------------------------------------------------
    # Phylogenetic logistic model
    # --------------------------------------------------------

    cat("Running phyloglm logistic_MPLE...\n")

    fit <- tryCatch(
        phyloglm(
            y ~ 1,
            data = dat,
            phy = tree,
            method = "logistic_MPLE"
        ),
        error = function(e) NULL
    )

    if (is.null(fit)) {

        cat("STATUS: MODEL FAILED\n")

        results[[length(results) + 1]] <- data.frame(
            class = cls,
            present = present,
            absent = absent,
            n = length(y),
            estimate = NA_real_,
            SE = NA_real_,
            z = NA_real_,
            p = NA_real_,
            alpha = NA_real_,
            status = "failed",
            stringsAsFactors = FALSE
        )

        next
    }

    sm <- summary(fit)

    coef_table <- sm$coefficients

    estimate <- unname(coef_table[1, "Estimate"])
    SE       <- unname(coef_table[1, "StdErr"])
    z        <- unname(coef_table[1, "z.value"])
    p        <- unname(coef_table[1, "p.value"])

    alpha <- unname(fit$alpha)

    cat("Estimate:", estimate, "\n")
    cat("SE      :", SE, "\n")
    cat("z       :", z, "\n")
    cat("P-value :", p, "\n")
    cat("Alpha   :", alpha, "\n")

    # 95% CI
    lower <- estimate - 1.96 * SE
    upper <- estimate + 1.96 * SE

    # Convert intercept to probability
    probability <- plogis(estimate)

    results[[length(results) + 1]] <- data.frame(
        class = cls,
        present = present,
        absent = absent,
        n = length(y),
        prevalence = present / length(y),
        estimate = estimate,
        SE = SE,
        CI_lower = lower,
        CI_upper = upper,
        probability = probability,
        z = z,
        p = p,
        alpha = alpha,
        status = "fitted",
        stringsAsFactors = FALSE
    )
}

# ============================================================
# COMBINE RESULTS
# ============================================================

res <- do.call(rbind, results)

# ============================================================
# BH-FDR
# ============================================================

res$FDR <- NA_real_

eligible <- res$status == "fitted" & !is.na(res$p)

res$FDR[eligible] <- p.adjust(
    res$p[eligible],
    method = "BH"
)

res$significant_FDR <- FALSE

res$significant_FDR[
    eligible & res$FDR < 0.05
] <- TRUE

# ============================================================
# SORT
# ============================================================

res <- res[order(
    res$FDR,
    res$p,
    na.last = TRUE
), ]

rownames(res) <- NULL

# ============================================================
# OUTPUT
# ============================================================

outfile <- "04B_phylogenetic_signal_main.tsv"

write.table(
    res,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("PRIMARY PHYLOGENETIC SIGNAL RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

print(
    res[
        res$status == "fitted",
        c(
            "class",
            "present",
            "absent",
            "prevalence",
            "estimate",
            "CI_lower",
            "CI_upper",
            "p",
            "FDR",
            "alpha",
            "significant_FDR"
        )
    ],
    row.names = FALSE
)

cat("\n")
cat("Eligible classes:", sum(res$status == "fitted"), "\n")
cat("Excluded classes:", sum(res$status == "excluded"), "\n")
cat("Failed classes  :", sum(res$status == "failed"), "\n")
cat(
    "FDR-significant:",
    sum(res$significant_FDR, na.rm = TRUE),
    "\n"
)

cat("\nOutput:\n")
cat(normalizePath(outfile), "\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("4B COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
