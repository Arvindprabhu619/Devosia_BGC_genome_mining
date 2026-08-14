## ============================================================
## Devosia genome metadata curation
##
## Purpose:
##   Curate NCBI Devosia assembly metadata by:
##   1. merging GCA/GCF paired assemblies
##   2. prioritizing RefSeq assemblies
##   3. excluding MAG/uncultured/binned assemblies
##   4. excluding highly fragmented assemblies (>200 scaffolds)
##
## Input:
##   data/metadata/Devosia.tsv
##
## Outputs:
##   data/metadata/devosia_curated_full.tsv
##   data/metadata/devosia_curated_keep_only.tsv
##   data/metadata/devosia_accessions.txt
##
## The project root is inferred from the location of this script,
## so the workflow does not depend on a user-specific filesystem path.
## ============================================================

## ---- Determine repository root ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(
    sub("^--file=", "", file_arg[1]),
    mustWork = TRUE
  )
  project_root <- normalizePath(
    file.path(dirname(script_path), "..", ".."),
    mustWork = TRUE
  )
} else {
  project_root <- normalizePath(".", mustWork = TRUE)
}

metadata_dir <- file.path(project_root, "data", "metadata")
input_file <- file.path(metadata_dir, "Devosia.tsv")

if (!file.exists(input_file)) {
  stop(
    "Input metadata file not found:\n",
    input_file,
    "\n\nPlace the NCBI metadata file at data/metadata/Devosia.tsv"
  )
}

## ---- Read metadata ----
df <- read.delim(
  input_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = TRUE,
  quote = ""
)

expected_header <- c(
  "Assembly.Name",
  "Assembly.Accession",
  "Assembly.Paired.Assembly.Accession",
  "Organism.Name",
  "Organism.Infraspecific.Names.Breed",
  "Organism.Infraspecific.Names.Strain",
  "Organism.Infraspecific.Names.Cultivar",
  "Organism.Infraspecific.Names.Ecotype",
  "Organism.Infraspecific.Names.Isolate",
  "Organism.Infraspecific.Names.Sex",
  "Annotation.Name",
  "Assembly.Level",
  "Assembly.Release.Date",
  "WGS.project.accession",
  "Assembly.Stats.Number.of.Scaffolds"
)

stopifnot(
  "Header layout changed -- check column order before renaming!" =
    identical(names(df), expected_header)
)

names(df) <- c(
  "AssemblyName",
  "AssemblyAccession",
  "PairedAccession",
  "OrganismName",
  "Breed",
  "Strain",
  "Cultivar",
  "Ecotype",
  "Isolate",
  "Sex",
  "AnnotationName",
  "AssemblyLevel",
  "ReleaseDate",
  "WGSAccession",
  "NumScaffolds"
)

cat("Total input rows:", nrow(df), "\n")

## ---- Merge GCA/GCF pairs ----
make_key <- function(acc, paired) {
  paired <- trimws(paired)

  if (nchar(paired) == 0 || is.na(paired)) {
    return(acc)
  }

  paste(sort(c(acc, paired)), collapse = "|")
}

df$GroupKey <- mapply(
  make_key,
  df$AssemblyAccession,
  df$PairedAccession
)

df$IsRefSeq <- grepl("^GCF_", df$AssemblyAccession)

groups <- split(df, df$GroupKey)

curated_rows <- lapply(groups, function(g) {

  has_refseq <- any(g$IsRefSeq)

  chosen <- if (has_refseq) {
    g[g$IsRefSeq, ][1, ]
  } else {
    g[1, ]
  }

  strain_lbl <- if (
    nchar(trimws(chosen$Strain)) > 0
  ) {
    trimws(chosen$Strain)
  } else {
    trimws(chosen$Isolate)
  }

  n_scaf <- suppressWarnings(
    as.integer(trimws(chosen$NumScaffolds))
  )

  data.frame(
    Assembly_Accession = chosen$AssemblyAccession,
    Has_RefSeq = ifelse(has_refseq, "yes", "no"),
    Organism = trimws(chosen$OrganismName),
    Strain = strain_lbl,
    Assembly_Level = trimws(chosen$AssemblyLevel),
    N_Scaffolds = ifelse(
      is.na(n_scaf),
      trimws(chosen$NumScaffolds),
      n_scaf
    ),
    Release_Date = trimws(chosen$ReleaseDate),
    AnnotationName = trimws(chosen$AnnotationName),
    stringsAsFactors = FALSE
  )
})

curated <- do.call(rbind, curated_rows)
rownames(curated) <- NULL

cat(
  "Unique physical assemblies (after GCA/GCF merge):",
  nrow(curated),
  "\n"
)

n_singleton_no_refseq <- sum(curated$Has_RefSeq == "no")

cat(
  "Assemblies with NO RefSeq pair in this dataset:",
  n_singleton_no_refseq,
  "\n\n"
)

## ---- Flag MAGs / uncultured / Candidatus / fragmented ----
mag_pattern <- paste0(
  "(?i)(MAG|bin\\.?\\d|binner|_bin|bin_|",
  "metabat|concoct|maxbin|refined|binchicken)"
)

is_mag_like <-
  grepl(mag_pattern, curated$Strain, perl = TRUE) |
  grepl(mag_pattern, curated$AnnotationName, perl = TRUE)

is_uncultured_or_candidatus <-
  grepl(
    "^(uncultured|candidatus)",
    curated$Organism,
    ignore.case = TRUE
  )

n_scaf_numeric <- suppressWarnings(
  as.integer(curated$N_Scaffolds)
)

highly_fragmented <-
  !is.na(n_scaf_numeric) &
  n_scaf_numeric > 200

exclude <-
  is_mag_like |
  is_uncultured_or_candidatus |
  highly_fragmented

exclude_reason <- ifelse(
  is_mag_like | is_uncultured_or_candidatus,
  "MAG/uncultured/binned",
  ifelse(
    highly_fragmented,
    paste0(">200 scaffolds (", n_scaf_numeric, ")"),
    ""
  )
)

exclude_reason <- ifelse(
  (is_mag_like | is_uncultured_or_candidatus) &
    highly_fragmented,
  paste0(
    exclude_reason,
    "; >200 scaffolds (",
    n_scaf_numeric,
    ")"
  ),
  exclude_reason
)

curated$Decision <- ifelse(
  exclude,
  "Remove",
  "Keep"
)

curated$Exclude_Reason <- exclude_reason

## ---- Confirm previously missed MAG ----
flagged <- curated$Decision[
  curated$Assembly_Accession == "GCA_056095215.1"
]

cat(
  "GCA_056095215.1 decision (should be 'Remove'):",
  flagged,
  "\n\n"
)

## ---- Rank + sort ----
level_rank <- c(
  "Complete Genome" = 0,
  "Chromosome" = 1,
  "Scaffold" = 2,
  "Contig" = 3
)

curated$Level_Rank <- ifelse(
  curated$Assembly_Level %in% names(level_rank),
  level_rank[curated$Assembly_Level],
  9
)

curated <- curated[
  order(
    curated$Decision != "Keep",
    curated$Organism,
    curated$Level_Rank,
    n_scaf_numeric
  ),
]

## ---- Write outputs ----
out_cols <- c(
  "Assembly_Accession",
  "Has_RefSeq",
  "Organism",
  "Strain",
  "Assembly_Level",
  "N_Scaffolds",
  "Release_Date",
  "Decision",
  "Exclude_Reason"
)

full_output <- file.path(
  metadata_dir,
  "devosia_curated_full.tsv"
)

keep_output <- file.path(
  metadata_dir,
  "devosia_curated_keep_only.tsv"
)

accession_output <- file.path(
  metadata_dir,
  "devosia_accessions.txt"
)

write.table(
  curated[, out_cols],
  full_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

kept <- curated[
  curated$Decision == "Keep",
  out_cols
]

write.table(
  kept,
  keep_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  kept$Assembly_Accession,
  accession_output
)

## ---- Summary ----
cat("=== SUMMARY ===\n")
cat("Kept:  ", nrow(kept), "\n")
cat(
  "Removed:",
  sum(curated$Decision == "Remove"),
  "\n"
)

cat(
  "Unique species among kept genomes:",
  length(unique(kept$Organism)),
  "\n\n"
)

## ---- Removal category breakdown ----
category <- ifelse(
  curated$Decision != "Remove",
  "Kept",
  ifelse(
    grepl(
      "MAG/uncultured/binned",
      curated$Exclude_Reason
    ) &
      grepl(
        ">200 scaffolds",
        curated$Exclude_Reason
      ),
    "MAG/uncultured + fragmented",
    ifelse(
      grepl(
        "MAG/uncultured/binned",
        curated$Exclude_Reason
      ),
      "MAG/uncultured only",
      "Fragmented only (>200 scaffolds)"
    )
  )
)

cat("=== Removal category breakdown ===\n")

print(
  table(
    category[curated$Decision == "Remove"]
  )
)

cat(
  "\nFiles written to:",
  metadata_dir,
  "\n"
)
