# Devosia BGC Genome Mining

**Genome-scale comparative genome mining of the genus *Devosia*: biosynthetic diversity, ecological association, and evolutionary stratification of secondary metabolism.**

---

## Overview

This repository contains the complete, reproducible computational pipeline for characterizing the biosynthetic gene cluster (BGC) repertoire of the bacterial genus *Devosia* from publicly available genomes.

**Key features:**

- 15-stage pipeline covering genome retrieval, quality control, taxonomy validation, phylogeny, BGC prediction, GCF analysis, statistics
- One-command reproducibility: bash run_pipeline.sh
- Publication-grade figures generated with unified palette
- Complete documentation in docs/

**Findings (from the manuscript):**

- Taxonomically validated Devosia/Devosia_A genome set
- BGCs across multiple biosynthetic classes
- ~82% putatively novel (no detectable MIBiG match)
- Rare GCFs dominate (89-96% across cutoffs)
- Niche-specific and evolutionarily stratified BGCs

---

## Quick Start

### Prerequisites

- Linux HPC with SLURM or equivalent
- 16+ cores, 128+ GB RAM
- 500 GB free disk space
- micromamba or conda

### Installation

    git clone https://github.com/Arvindprabhu619/Devosia_BGC_genome_mining.git
    cd Devosia_BGC_genome_mining
    micromamba create -f environment.yml
    micromamba activate devosia_bgc
    bash scripts/utils/download_databases.sh
    bash scripts/utils/check_deps.sh

### Run the full pipeline

    bash run_pipeline.sh              # all 15 stages
    bash run_pipeline.sh 06 10        # stages 6-10
    bash run_stage.sh 06              # just stage 6

---

## Pipeline Stages

| # | Stage | Purpose | Runtime |
|---|-------|---------|---------|
| 01 | Download genomes | Fetch Devosia assemblies from NCBI | 15-30 min |
| 02 | Curate metadata | Remove MAGs, duplicates, fragmented | 2 min |
| 03 | Build manifest | Map accessions to files | 2 min |
| 04 | CheckM2 | Genome completeness/contamination | 20-60 min |
| 05 | GTDB-Tk classify | Taxonomic validation (GTDB R232) | 20-60 min |
| 05b | Build phylogeny | bac120 alignment + FastTree | 30-60 min |
| 06 | antiSMASH | BGC prediction | 20-24 hr |
| 07 | Extract BGCs | Parse antiSMASH JSON | 5 min |
| 08 | Classify novelty | KnownClusterBlast vs MIBiG | 5 min |
| 09 | BiG-SCAPE | GCF clustering (c0.3/0.5/0.7) | 2-4 hr |
| 10 | GCF prevalence | Family distribution across genomes | 15 min |
| 11 | Phylogenetic stats | Pagel lambda, logistic regression | 1-2 hr |
| 12 | Prioritize | Score BGCs for follow-up | 15 min |
| 13 | Make figures | Publication-grade figures | 30 min |
| 14 | Ecological | Niche association (Fisher exact) | 30 min |
| 15 | Phylostratigraphy | GCF evolutionary ages | 1 hr |

Total runtime: ~30-40 hours on a 16-core, 128 GB node.

---

## Terminology (Important)

- **Putatively novel BGC** — no detectable KnownClusterBlast match to MIBiG v3.1. Does NOT prove novel chemistry.
- **Related BGC** — detectable MIBiG match at <70% gene similarity.
- **Known BGC** — detectable MIBiG match at >=70% gene similarity. Study-specific cutoff.
- **GCF** — Gene Cluster Family
- **Rare GCF** — prevalence <10% of genomes
- **Core GCF** — prevalence >=90% of genomes

---

## Tools and Versions

| Tool | Version |
|------|---------|
| antiSMASH | 7.1.0 |
| CheckM2 | 1.1.0 |
| GTDB-Tk | 2.7.2 |
| BiG-SCAPE | 2.0.3 |
| FastTree | 2.2.0 |
| R | 4.4+ |
| MIBiG | 3.1 |
| GTDB | R232 |

---

## Reproducibility

- Every stage logs to logs/
- Output tables include SHA256 checksums
- Full environment captured in environment.yml
- Database versions in docs/04_databases.md

---

## Citation

Manuscript in preparation.

---

## Contact

- **Corresponding author:** Dr Inamul Hasan Madar (CIODS, Yenepoya University)
- **Pipeline maintainer:** Arvind Prabhu
- **Issues:** Open a GitHub issue

---

## License

MIT License. See LICENSE.
