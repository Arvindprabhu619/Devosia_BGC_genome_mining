###############################################################################
# ANALYSIS 4A
# PHYLOGENETIC LOGISTIC REGRESSION:
# ZERO-LENGTH BRANCH SENSITIVITY ANALYSIS
#
# Dataset:
#   110 Devosia genomes
#   Table S2: 21 BGC classes
#
# Sensitivity:
#   epsilon = 1e-8
#   epsilon = 1e-7
#   epsilon = 1e-6
#
# Model:
#   phyloglm(..., method = "logistic_MPLE")
#
# Multiple testing:
#   Benjamini-Hochberg FDR within each epsilon
#
# IMPORTANT:
#   Sensitivity analysis only.
###############################################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
})

cat("\n")
cat(strrep("=", 80), "\n")
cat("ANALYSIS 4A: PHYLOGLM ZERO-LENGTH BRANCH SENSITIVITY\n")
cat(strrep("=", 80), "\n")

cat("Working directory:", getwd(), "\n")
cat("R version:", R.version.string, "\n")
cat("phylolm version:", as.character(packageVersion("phylolm")), "\n")


###############################################################################
# INPUTS
###############################################################################

tree_file <- "../phylogeny_110/Devosia_110_final_rooted_pruned.tree"

s2_file <- "../manuscript_outputs/TableS2_genome_by_class_matrix.tsv"

output_file <- "04A_phyloglm_zero_branch_sensitivity.tsv"


###############################################################################
# READ TREE
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("READING ORIGINAL TREE\n")
cat(strrep("=", 80), "\n")

tree_original <- read.tree(tree_file)

cat("Tree tips:", length(tree_original$tip.label), "\n")
cat("Rooted:", is.rooted(tree_original), "\n")
cat("Branches:", length(tree_original$edge.length), "\n")
cat("Zero-length branches:",
    sum(tree_original$edge.length <= 0), "\n")
cat("Minimum branch length:",
    min(tree_original$edge.length), "\n")
cat("Maximum branch length:",
    max(tree_original$edge.length), "\n")


###############################################################################
# READ TABLE S2
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("READING TABLE S2\n")
cat(strrep("=", 80), "\n")

s2 <- read.delim(
  s2_file,
  sep = "\t",
  check.names = FALSE
)

cat("Rows:", nrow(s2), "\n")
cat("Columns:", ncol(s2), "\n")


###############################################################################
# CANONICAL ACCESSION EXTRACTION
###############################################################################

# Tree labels contain:
#
# GCF_XXXXXXXXX.X_ASSEMBLY_NAME_genomic
#
# Extract only:
#
# GCF_XXXXXXXXX.X

tree_acc <- sub(
  "^((GC[FA]_[0-9]+\\.[0-9]+)).*$",
  "\\1",
  tree_original$tip.label
)

###############################################################################
# ACCESSION VALIDATION
###############################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("ACCESSION RECONCILIATION\n")
cat(strrep("-", 80), "\n")

cat("Tree:", length(tree_acc), "\n")
cat("Table S2:", nrow(s2), "\n")
cat("Common:", length(intersect(tree_acc, s2$accession)), "\n")

stopifnot(length(tree_acc) == 110)
stopifnot(length(unique(tree_acc)) == 110)

stopifnot(
  length(intersect(tree_acc, s2$accession)) == 110
)

stopifnot(
  length(setdiff(tree_acc, s2$accession)) == 0
)

stopifnot(
  length(setdiff(s2$accession, tree_acc)) == 0
)

cat("PASS: All 110 genomes match.\n")


###############################################################################
# ALIGN TABLE S2 TO TREE
###############################################################################

s2 <- s2[match(tree_acc, s2$accession), ]

stopifnot(all(s2$accession == tree_acc))

cat("PASS: Table S2 aligned to tree order.\n")


###############################################################################
# BGC CLASSES
###############################################################################

bgc_classes <- setdiff(colnames(s2), "accession")

cat("\n")
cat(strrep("=", 80), "\n")
cat("BGC CLASSES\n")
cat(strrep("=", 80), "\n")

cat("Total classes:", length(bgc_classes), "\n")
print(bgc_classes)


###############################################################################
# SENSITIVITY EPSILONS
###############################################################################

epsilons <- c(
  1e-8,
  1e-7,
  1e-6
)


###############################################################################
# PREVALENCE THRESHOLD
###############################################################################

# Classes with <5 or >105 genomes are excluded because the binary response
# is too sparse / constant for this analysis.

min_present <- 5
min_absent <- 5


###############################################################################
# STORAGE
###############################################################################

results <- list()

counter <- 0


###############################################################################
# LOOP OVER EPSILONS
###############################################################################

for (eps in epsilons) {
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("EPSILON =", format(eps, scientific = TRUE), "\n")
  cat(strrep("=", 80), "\n")
  
  tree <- tree_original
  
  zero_before <- sum(tree$edge.length <= 0)
  
  tree$edge.length[tree$edge.length <= 0] <- eps
  
  zero_after <- sum(tree$edge.length <= 0)
  
  min_branch <- min(tree$edge.length)
  
  cat("Zero branches before:", zero_before, "\n")
  cat("Zero branches after :", zero_after, "\n")
  cat("Minimum branch      :", min_branch, "\n")
  
  
  ###########################################################################
  # LOOP OVER BGC CLASSES
  ###########################################################################
  
  epsilon_results <- list()
  
  for (class_name in bgc_classes) {
    
    cat("\nCLASS:", class_name, "\n")
    
    # Binary phenotype
    y <- as.numeric(s2[[class_name]] > 0)
    
    n_present <- sum(y == 1, na.rm = TRUE)
    n_absent  <- sum(y == 0, na.rm = TRUE)
    n_na      <- sum(is.na(y))
    
    prevalence <- 100 * n_present / length(y)
    
    cat("Present:", n_present, "\n")
    cat("Absent :", n_absent, "\n")
    
    #######################################################################
    # PREVALENCE FILTER
    #######################################################################
    
    if (
      n_na > 0 ||
      n_present < min_present ||
      n_absent < min_absent
    ) {
      
      cat("STATUS: EXCLUDED due to prevalence/NA\n")
      
      counter <- counter + 1
      
      epsilon_results[[length(epsilon_results) + 1]] <- data.frame(
        epsilon = eps,
        class = class_name,
        genomes_total = length(y),
        genomes_with_class = n_present,
        genomes_without_class = n_absent,
        prevalence_percent = prevalence,
        status = "excluded",
        estimate = NA_real_,
        SE = NA_real_,
        z = NA_real_,
        P_value = NA_real_,
        FDR_BH = NA_real_,
        alpha = NA_real_,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    
    #######################################################################
    # CRITICAL DATA/TREE MATCH
    #######################################################################
    
    # Use EXACT tree tip labels as row names.
    #
    # This is essential for phylolm().
    #
    dat <- data.frame(
      y = y,
      row.names = tree$tip.label
    )
    
    stopifnot(
      identical(rownames(dat), tree$tip.label)
    )
    
    
    #######################################################################
    # FIT PHYLOGLM
    #######################################################################
    
    cat("Running phyloglm logistic_MPLE...\n")
    
    fit <- tryCatch(
      
      phyloglm(
        y ~ 1,
        data = dat,
        phy = tree,
        method = "logistic_MPLE"
      ),
      
      error = function(e) e
    )
    
    
    #######################################################################
    # HANDLE FAILURE
    #######################################################################
    
    if (inherits(fit, "error")) {
      
      cat("FAILED:", conditionMessage(fit), "\n")
      
      counter <- counter + 1
      
      epsilon_results[[length(epsilon_results) + 1]] <- data.frame(
        epsilon = eps,
        class = class_name,
        genomes_total = length(y),
        genomes_with_class = n_present,
        genomes_without_class = n_absent,
        prevalence_percent = prevalence,
        status = "failed",
        estimate = NA_real_,
        SE = NA_real_,
        z = NA_real_,
        P_value = NA_real_,
        FDR_BH = NA_real_,
        alpha = NA_real_,
        stringsAsFactors = FALSE
      )
      
      next
    }
    
    
    #######################################################################
    # EXTRACT COEFFICIENT STATISTICS
    #######################################################################
    
    coef_matrix <- summary(fit)$coefficients
    
    estimate <- unname(coef_matrix[1, "Estimate"])
    
    SE <- unname(coef_matrix[1, "StdErr"])
    
    z <- unname(coef_matrix[1, "z.value"])
    
    P <- unname(coef_matrix[1, "p.value"])
    
    alpha <- unname(fit$alpha)
    
    
    #######################################################################
    # OUTPUT
    #######################################################################
    
    cat("Estimate:", estimate, "\n")
    cat("SE      :", SE, "\n")
    cat("z       :", z, "\n")
    cat("P-value :", P, "\n")
    cat("Alpha   :", alpha, "\n")
    
    counter <- counter + 1
    
    epsilon_results[[length(epsilon_results) + 1]] <- data.frame(
      epsilon = eps,
      class = class_name,
      genomes_total = length(y),
      genomes_with_class = n_present,
      genomes_without_class = n_absent,
      prevalence_percent = prevalence,
      status = "successful",
      estimate = estimate,
      SE = SE,
      z = z,
      P_value = P,
      FDR_BH = NA_real_,
      alpha = alpha,
      stringsAsFactors = FALSE
    )
  }
  
  
  ###########################################################################
  # BH-FDR WITHIN EPSILON
  ###########################################################################
  
  epsilon_df <- do.call(
    rbind,
    epsilon_results
  )
  
  valid <- (
    epsilon_df$status == "successful" &
      !is.na(epsilon_df$P_value)
  )
  
  epsilon_df$FDR_BH[valid] <- p.adjust(
    epsilon_df$P_value[valid],
    method = "BH"
  )
  
  results[[length(results) + 1]] <- epsilon_df
}


###############################################################################
# COMBINE
###############################################################################

final_results <- do.call(
  rbind,
  results
)

rownames(final_results) <- NULL


###############################################################################
# WRITE OUTPUT
###############################################################################

write.table(
  final_results,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("\n")
cat(strrep("=", 80), "\n")
cat("OUTPUT WRITTEN\n")
cat(strrep("=", 80), "\n")
cat(output_file, "\n")


###############################################################################
# SUMMARY
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("SUMMARY BY EPSILON\n")
cat(strrep("=", 80), "\n")

for (eps in epsilons) {
  
  x <- final_results[
    final_results$epsilon == eps,
  ]
  
  tested <- sum(x$status == "successful")
  
  failed <- sum(x$status == "failed")
  
  significant <- sum(
    x$status == "successful" &
      !is.na(x$FDR_BH) &
      x$FDR_BH < 0.05
  )
  
  cat("\n")
  cat("epsilon =", format(eps, scientific = TRUE), "\n")
  cat("Tested:", tested, "\n")
  cat("Failed:", failed, "\n")
  cat("Significant after BH-FDR:", significant, "\n")
}


###############################################################################
# SENSITIVITY COMPARISON
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("SENSITIVITY COMPARISON\n")
cat(strrep("=", 80), "\n")

successful <- final_results[
  final_results$status == "successful",
]

if (nrow(successful) > 0) {
  
  for (class_name in unique(successful$class)) {
    
    x <- successful[
      successful$class == class_name,
    ]
    
    if (nrow(x) == length(epsilons)) {
      
      est_range <- max(x$estimate) - min(x$estimate)
      
      p_range <- max(x$P_value) - min(x$P_value)
      
      fdr_range <- max(x$FDR_BH) - min(x$FDR_BH)
      
      cat(
        sprintf(
          "%-25s estimate range = %.6g; P range = %.6g; FDR range = %.6g\n",
          class_name,
          est_range,
          p_range,
          fdr_range
        )
      )
    }
  }
}


###############################################################################
# FINAL STATUS
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("FINAL STATUS\n")
cat(strrep("=", 80), "\n")

total_success <- sum(
  final_results$status == "successful"
)

total_failed <- sum(
  final_results$status == "failed"
)

cat("Total model fits attempted:", total_success + total_failed, "\n")
cat("Successful model fits:", total_success, "\n")
cat("Failed model fits:", total_failed, "\n")

if (total_failed == 0) {
  cat("\nPASS\n")
  cat("All eligible binary BGC classes successfully fitted\n")
  cat("across all epsilon values.\n")
} else {
  cat("\nWARNING\n")
  cat("Some model fits failed.\n")
}

cat("\nOutput:\n")
cat("  ", output_file, "\n")

cat(strrep("=", 80), "\n")
RSCRIPT
