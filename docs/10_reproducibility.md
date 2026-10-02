# Reproducibility

## Verification protocol

### 1. Environment verification

    bash scripts/utils/check_deps.sh

Checks all tool wrappers are functional and databases exist.

### 2. Database checksums

SHA256 checksums recorded in env/database_versions.txt.

    sha256sum -c env/database_versions.txt

### 3. Stage-by-stage verification

Each stage writes checksums to results/<stage>/checksums.sha256:

    cd results/06_antismash/parsed
    sha256sum -c checksums.sha256

### 4. End-to-end reproducibility

For full reproducibility, run the pipeline from a clean state:

    rm -rf results/ figures/ tables/
    bash run_pipeline.sh

Compare outputs against published results.

## Determinism

The pipeline is designed to be deterministic:

- Same inputs -> same outputs
- Random seeds fixed where applicable
- Tool versions pinned in environment.yml
- Database versions documented in docs/04_databases.md

Note: antiSMASH and GTDB-Tk may produce slightly different results on different runs due to parallel processing and database updates. Pin database versions for full reproducibility.

## Software versions

    bash scripts/utils/check_deps.sh > env/software_versions.txt

Commit env/software_versions.txt to git to record the environment used for a published analysis.

## Logs

Every stage writes logs to logs/<stage>_<timestamp>.log.

Include logs when reporting issues.

## Reporting issues

Open a GitHub issue with:
- The failing stage number
- The command used
- The log file contents
- Your environment details (bash scripts/utils/check_deps.sh output)
