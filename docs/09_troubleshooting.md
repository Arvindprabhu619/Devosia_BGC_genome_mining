# Troubleshooting

## Common errors and fixes

### Broken base Python

Symptom: Fatal Python error: Failed to import encodings module

Cause: Broken conda base installation.

Fix: Use env-specific Python via $PYTHON_BIN (set to checkm2 env Python). All Python calls in pipeline scripts use $PYTHON_BIN.

### GTDB-Tk: reference data does not exist

Symptom: The GTDB-Tk reference data does not exist or is corrupted.

Cause: GTDBTK_DATA_PATH points to a non-existent path.

Fix: Export GTDBTK_DATA_PATH=$(GTDBTK_DATA). This is handled by scripts/00_paths.sh.

### GTDB-Tk: same id as GTDB-Tk reference genomes

Symptom: You have N genomes with the same id as GTDB-Tk reference genomes, please rename them.

Cause: Genome accessions match GTDB reference accessions.

Fix: Symlinks use query_ prefix. Handled automatically by scripts/05_run_gtdbtk.sh.

### GTDB-Tk align: FileNotFoundError marker_genes

Symptom: FileNotFoundError: ... intermediate_results/marker_genes

Cause: Passing wrong path to --identify_dir.

Fix: Use --identify_dir results/05b_phylogeny/identify_all/identify (not the parent dir).

### antiSMASH: unrecognized arguments --smcogs

Symptom: antismash: error: unrecognized arguments: --smcogs

Cause: --smcogs was removed in antiSMASH 7.x.

Fix: Do not use --smcogs. Removed from scripts/06_run_antismash.sh.

### antiSMASH: Output directory contains other files

Symptom: ERROR: Output directory contains other files, aborting for safety

Cause: antiSMASH refuses to write into a pre-populated directory.

Fix: Log files go to results/06_antismash/logs/ (separate from output dir). Handled by scripts/06_run_antismash.sh.

### BiG-SCAPE: no output

Symptom: BiG-SCAPE finishes but no cluster files

Cause: Pfam index not built, or input GBK files missing.

Fix: Ensure Pfam-A.hmm.h3f/.h3i/.h3m/.h3p exist (run hmmpress Pfam-A.hmm). Verify results/06_antismash/raw_output/*/*.gbk present.

### R: package not found

Symptom: Error in library(X) : there is no package called X

Fix: Install via micromamba install r-X (see environment.yml).

## Resume from a failed stage

    bash run_pipeline.sh <stage> 15

All scripts are idempotent - they skip completed work.
