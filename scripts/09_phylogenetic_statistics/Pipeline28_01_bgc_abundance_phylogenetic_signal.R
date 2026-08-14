###############################################################################
# ANALYSIS 1
# PHYLOGENETIC SIGNAL OF TOTAL BGC ABUNDANCE
#
# Dataset:
#   110 Devosia genomes
#   582 BGC regions
#
# Primary test:
#   Pagel's lambda
#
# R version:
#   Compatible with R 4.1.2
###############################################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

cat("\n")
cat(strrep("=", 80), "\n")
cat("ANALYSIS 1: PHYLOGENETIC SIGNAL OF TOTAL BGC ABUNDANCE\n")
cat(strrep("=", 80), "\n")

cat("Working directory:", getwd(), "\n")
cat("R version:", R.version.string, "\n")


###############################################################################
# FILES
###############################################################################

tree_file <- "../phylogeny_110/Devosia_110_final_rooted_pruned.tree"
bgc_file  <- "../bgc_counts_per_genome_fixed.tsv"

aligned_file <- "01_BGC_abundance_trait_aligned.tsv"
result_file  <- "01_BGC_abundance_Pagels_lambda.tsv"


###############################################################################
# ROBUST ACCESSION EXTRACTION
###############################################################################

canonical_accession <- function(x) {
  
  x <- as.character(x)
  
  # Initialize output
  out <- rep(NA_character_, length(x))
  
  # Identify strings containing GCF/GCA accession
  hit <- grepl(
    "GC[FA]_[0-9]+\\.[0-9]+",
    x,
    perl = TRUE
  )
  
  # Extract accession using sub()
  # This is robust in older R versions.
  out[hit] <- sub(
    ".*(GC[FA]_[0-9]+\\.[0-9]+).*",
    "\\1",
    x[hit],
    perl = TRUE
  )
  
  out
}


###############################################################################
# INPUT VALIDATION
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("INPUT FILE VALIDATION\n")
cat(strrep("=", 80), "\n")

cat("Tree file:", tree_file, "\n")
cat("BGC file :", bgc_file, "\n")

if (!file.exists(tree_file)) {
  stop("ERROR: Tree file does not exist: ", tree_file)
}

if (!file.exists(bgc_file)) {
  stop("ERROR: BGC file does not exist: ", bgc_file)
}


###############################################################################
# READ TREE
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("READING TREE\n")
cat(strrep("=", 80), "\n")

tree <- read.tree(tree_file)

cat("Tree tips:", Ntip(tree), "\n")
cat("Tree rooted:", is.rooted(tree), "\n")
cat("Number of branches:", Nedge(tree), "\n")

if (is.null(tree$edge.length)) {
  stop("ERROR: Tree has no branch lengths.")
}

cat("Branch lengths present: TRUE\n")

cat(
  "Minimum branch length:",
  min(tree$edge.length),
  "\n"
)

cat(
  "Maximum branch length:",
  max(tree$edge.length),
  "\n"
)

cat(
  "Zero-length branches:",
  sum(tree$edge.length == 0),
  "\n"
)

cat(
  "Negative branch lengths:",
  sum(tree$edge.length < 0),
  "\n"
)

if (any(tree$edge.length < 0)) {
  stop("ERROR: Tree contains negative branch lengths.")
}


###############################################################################
# CANONICALIZE TREE ACCESSIONS
###############################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("CANONICALIZING TREE ACCESSIONS\n")
cat(strrep("-", 80), "\n")

original_tree_labels <- tree$tip.label

tree_ids <- canonical_accession(original_tree_labels)

cat(
  "Original tree tip labels:",
  length(original_tree_labels),
  "\n"
)

cat(
  "Canonical tree accessions:",
  sum(!is.na(tree_ids)),
  "\n"
)

cat(
  "Unique canonical accessions:",
  length(unique(tree_ids)),
  "\n"
)

# Show failures if any
if (any(is.na(tree_ids))) {
  
  cat("\nERROR: Could not canonicalize these tree labels:\n")
  
  print(original_tree_labels[is.na(tree_ids)])
  
  stop("Tree accession canonicalization failed.")
}

# Check duplicates
if (anyDuplicated(tree_ids)) {
  
  cat("\nERROR: Duplicate canonical tree accessions detected:\n")
  
  print(unique(tree_ids[duplicated(tree_ids)]))
  
  stop("Tree contains duplicate accession IDs.")
}

if (length(tree_ids) != 110) {
  
  stop(
    "ERROR: Expected 110 tree genomes but found ",
    length(tree_ids),
    "."
  )
}

###############################################################################
# CRITICAL STEP
#
# Replace the tree's original labels with canonical accession IDs.
###############################################################################

tree$tip.label <- tree_ids

cat(
  "PASS: Tree contains",
  length(tree$tip.label),
  "unique canonical accession IDs.\n"
)


###############################################################################
# READ BGC TABLE
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("READING BGC TABLE\n")
cat(strrep("=", 80), "\n")

bgc <- read.table(
  bgc_file,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("accession", "BGC_count")
)

cat("Rows:", nrow(bgc), "\n")
cat("Columns:", ncol(bgc), "\n")


###############################################################################
# BGC ACCESSION VALIDATION
###############################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("BGC TABLE VALIDATION\n")
cat(strrep("-", 80), "\n")

bgc$accession_original <- bgc$accession

bgc$accession <- canonical_accession(
  bgc$accession
)

missing_accessions <- sum(
  is.na(bgc$accession)
)

missing_counts <- sum(
  is.na(bgc$BGC_count)
)

duplicate_accessions <- sum(
  duplicated(bgc$accession)
)

cat(
  "Missing accessions:",
  missing_accessions,
  "\n"
)

cat(
  "Missing BGC counts:",
  missing_counts,
  "\n"
)

cat(
  "Duplicate accessions:",
  duplicate_accessions,
  "\n"
)

if (missing_accessions > 0) {
  stop("ERROR: Invalid BGC accession detected.")
}

if (missing_counts > 0) {
  stop("ERROR: Missing BGC count detected.")
}

if (duplicate_accessions > 0) {
  stop("ERROR: Duplicate BGC accession detected.")
}

if (nrow(bgc) != 110) {
  stop(
    "ERROR: Expected 110 BGC genomes but found ",
    nrow(bgc),
    "."
  )
}

cat("PASS: BGC table contains 110 unique valid accessions.\n")


###############################################################################
# BGC COUNT VALIDATION
###############################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("BGC COUNT VALIDATION\n")
cat(strrep("-", 80), "\n")

bgc$BGC_count <- as.numeric(bgc$BGC_count)

cat(
  "Minimum BGC count:",
  min(bgc$BGC_count),
  "\n"
)

cat(
  "Maximum BGC count:",
  max(bgc$BGC_count),
  "\n"
)

cat(
  "Mean BGC count:",
  mean(bgc$BGC_count),
  "\n"
)

cat(
  "Median BGC count:",
  median(bgc$BGC_count),
  "\n"
)

cat(
  "Total BGCs:",
  sum(bgc$BGC_count),
  "\n"
)

if (any(bgc$BGC_count < 0)) {
  stop("ERROR: Negative BGC counts detected.")
}

if (sum(bgc$BGC_count) != 582) {
  stop(
    "ERROR: Expected 582 total BGCs but found ",
    sum(bgc$BGC_count),
    "."
  )
}

cat("PASS: Total BGC count = 582.\n")


###############################################################################
# TREE ↔ BGC ACCESSION RECONCILIATION
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("TREE ↔ BGC ACCESSION RECONCILIATION\n")
cat(strrep("=", 80), "\n")

tree_set <- unique(tree$tip.label)
bgc_set  <- unique(bgc$accession)

common_ids <- intersect(
  tree_set,
  bgc_set
)

tree_only <- setdiff(
  tree_set,
  bgc_set
)

bgc_only <- setdiff(
  bgc_set,
  tree_set
)

cat("Tree genomes:", length(tree_set), "\n")
cat("BGC genomes:", length(bgc_set), "\n")
cat("Common:", length(common_ids), "\n")
cat("Tree only:", length(tree_only), "\n")
cat("BGC only:", length(bgc_only), "\n")

if (length(tree_only) > 0) {
  
  cat("\nTree-only accessions:\n")
  print(tree_only)
  
  stop("ERROR: Tree/BGC mismatch.")
}

if (length(bgc_only) > 0) {
  
  cat("\nBGC-only accessions:\n")
  print(bgc_only)
  
  stop("ERROR: Tree/BGC mismatch.")
}

if (length(common_ids) != 110) {
  stop("ERROR: Expected 110 common genomes.")
}

cat("PASS: All 110 accessions match exactly.\n")


###############################################################################
# MATCH BGC TRAIT TO TREE
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("MATCHING BGC TRAIT TO TREE\n")
cat(strrep("=", 80), "\n")

trait <- bgc$BGC_count

names(trait) <- bgc$accession

# Exact tree order
trait <- trait[
  match(tree$tip.label, names(trait))
]

names(trait) <- tree$tip.label

cat(
  "Trait values matched:",
  length(trait),
  "\n"
)

cat(
  "Trait names:",
  length(names(trait)),
  "\n"
)

cat(
  "Missing trait values:",
  sum(is.na(trait)),
  "\n"
)

if (length(trait) != Ntip(tree)) {
  stop("ERROR: Trait length does not equal tree tips.")
}

if (any(is.na(trait))) {
  stop("ERROR: Missing trait values.")
}

if (!identical(names(trait), tree$tip.label)) {
  stop("ERROR: Trait names are not identical to tree tip labels.")
}

cat("PASS: Trait is correctly aligned to the phylogenetic tree.\n")


###############################################################################
# SAVE ALIGNED TRAIT
###############################################################################

aligned <- data.frame(
  accession = tree$tip.label,
  BGC_count = as.numeric(trait)
)

write.table(
  aligned,
  aligned_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Saved:",
  aligned_file,
  "\n"
)


###############################################################################
# FINAL TRAIT CHECK
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("FINAL TRAIT CHECK\n")
cat(strrep("=", 80), "\n")

cat(
  "Number of trait values:",
  length(trait),
  "\n"
)

cat(
  "Number of tree tips:",
  Ntip(tree),
  "\n"
)

cat(
  "Total BGCs represented:",
  sum(trait),
  "\n"
)

cat(
  "Minimum:",
  min(trait),
  "\n"
)

cat(
  "Maximum:",
  max(trait),
  "\n"
)

cat(
  "Mean:",
  mean(trait),
  "\n"
)

cat(
  "Median:",
  median(trait),
  "\n"
)

if (
  length(trait) != 110 ||
  Ntip(tree) != 110 ||
  sum(trait) != 582
) {
  stop("ERROR: Final trait validation failed.")
}

cat("PASS: Tree and BGC trait are correctly aligned.\n")


###############################################################################
# FINAL PHYLOGENETIC ALIGNMENT CHECK
###############################################################################

cat("\n")
cat(strrep("-", 80), "\n")
cat("PHYLOGENETIC ALIGNMENT CHECK\n")
cat(strrep("-", 80), "\n")

cat(
  "Tree tip 1:",
  tree$tip.label[1],
  "\n"
)

cat(
  "Trait name 1:",
  names(trait)[1],
  "\n"
)

cat(
  "Tree tip 110:",
  tree$tip.label[110],
  "\n"
)

cat(
  "Trait name 110:",
  names(trait)[110],
  "\n"
)

if (!all(tree$tip.label == names(trait))) {
  stop("ERROR: Tree and trait labels are not identical.")
}

cat(
  "PASS: All 110 tree labels exactly match trait names.\n"
)


###############################################################################
# PAGEL'S LAMBDA
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("PAGEL'S LAMBDA\n")
cat(strrep("=", 80), "\n")

cat(
  "Testing whether total BGC abundance shows phylogenetic signal.\n"
)

cat(
  "Null hypothesis: lambda = 0 (no phylogenetic signal).\n"
)

cat(
  "Alternative hypothesis: lambda > 0 (phylogenetic covariance is present).\n"
)

cat(
  "\nRunning phytools::phylosig()...\n"
)

lambda_result <- tryCatch(
  
  phylosig(
    tree,
    trait,
    method = "lambda",
    test = TRUE
  ),
  
  error = function(e) {
    
    cat("\nERROR during Pagel lambda calculation:\n")
    cat(conditionMessage(e), "\n")
    
    stop(
      "Pagel lambda calculation failed."
    )
  }
)


###############################################################################
# EXTRACT RESULTS
###############################################################################

lambda_estimate <- as.numeric(
  lambda_result$lambda
)

log_likelihood <- as.numeric(
  lambda_result$logL
)

p_value <- as.numeric(
  lambda_result$P
)

cat("\n")
cat(
  "Phylogenetic signal lambda :",
  format(lambda_estimate, digits = 10),
  "\n"
)

cat(
  "logL(lambda)               :",
  format(log_likelihood, digits = 10),
  "\n"
)

cat(
  "P-value                    :",
  format(p_value, digits = 10),
  "\n"
)


###############################################################################
# INTERPRETATION
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("INTERPRETATION\n")
cat(strrep("=", 80), "\n")

if (is.na(p_value)) {
  
  cat(
    "WARNING: P-value could not be calculated.\n"
  )
  
} else if (p_value < 0.05) {
  
  cat(
    "RESULT: Significant phylogenetic signal detected.\n\n"
  )
  
  cat(
    "Total BGC abundance is significantly associated ",
    "with phylogenetic relatedness among the 110 ",
    "Devosia genomes.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "RESULT: No significant phylogenetic signal detected.\n\n"
  )
  
  cat(
    "Total BGC abundance does not show statistically ",
    "significant association with phylogenetic ",
    "relatedness among the 110 Devosia genomes.\n",
    sep = ""
  )
}


###############################################################################
# SAVE RESULTS
###############################################################################

result <- data.frame(
  analysis = "Phylogenetic signal of total BGC abundance",
  trait = "Total_BGC_abundance",
  n_genomes = length(trait),
  total_BGCs = sum(trait),
  mean_BGCs = mean(trait),
  median_BGCs = median(trait),
  min_BGCs = min(trait),
  max_BGCs = max(trait),
  Pagel_lambda = lambda_estimate,
  logLik_lambda = log_likelihood,
  P_value = p_value,
  stringsAsFactors = FALSE
)

write.table(
  result,
  result_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")
cat(
  "Saved result:",
  result_file,
  "\n"
)


###############################################################################
# FINAL STATUS
###############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("FINAL STATUS\n")
cat(strrep("=", 80), "\n")

cat("PASS\n")
cat("✓ 110 tree genomes\n")
cat("✓ 110 BGC genomes\n")
cat("✓ Exact accession matching\n")
cat("✓ 582 total BGCs\n")
cat("✓ Trait aligned to tree\n")
cat("✓ Pagel's lambda calculated\n")
cat("✓ Results saved\n")

cat("\nOutput files:\n")
cat("  ", aligned_file, "\n", sep = "")
cat("  ", result_file, "\n", sep = "")

cat("\n")
cat(strrep("=", 80), "\n")

RSCRIPT
