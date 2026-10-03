# Devosia BGC Genome Mining

**Genome-scale comparative genome mining of the genus *Devosia*: biosynthetic diversity, ecological association, and evolutionary stratification of secondary metabolism.**

[![Pipeline](https://img.shields.io/badge/pipeline-15_stages-blue)]()
[![Reproducible](https://img.shields.io/badge/reproducible-yes-brightgreen)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

---

## Overview

This repository contains the complete, reproducible computational pipeline
for characterizing the biosynthetic gene cluster (BGC) repertoire of the
bacterial genus *Devosia* from publicly available genomes.

**Key features:**

- 15-stage pipeline covering genome retrieval, quality control, taxonomy validation, phylogeny, BGC prediction, GCF analysis, statistics
- One-command reproducibility: bash run_pipeline.sh
- Publication-grade figures generated with unified palette
- Complete documentation in docs/

---

## Key Results

| Metric | Value |
|--------|-------|
| Final genomes (Devosia + Devosia_A) | 124 |
| Physical BGC regions | 663 |
| Hybrid BGC regions | 27 |
| Biosynthetic class assignments | 690 |
| Distinct biosynthetic classes | 21 |
| Putatively novel BGCs | 538 (81.15%) |
| Related BGCs (less than 70% KCB) | 81 (12.22%) |
| Known BGCs (70% or more KCB) | 44 (6.64%) |
| GCFs at c0.3 / c0.5 / c0.7 | 217 / 145 / 113 |
| Core GCFs (90% or more prevalence) | 1 |
| Rare GCFs (less than 10% prevalence) | 89-96% across cutoffs |
| Devosia_A-specific GCFs | 26 |
| Shared Devosia/Devosia_A GCFs | 14 |
| Pagel's lambda (BGC abundance) | 0.933 (P = 3.2e-11) |
| Pagel's lambda (novelty proportion) | 0.981 (P = 2.5e-29) |
| Significant niche-BGC associations (FDR<0.05) | 5 |
| Ancient GCFs | 46 |
| Singleton GCFs | 38 |
| Prioritized candidate BGCs | 7 |

**Highlights:**

- Devosia genomes encode 663 BGCs across 21 biosynthetic classes
- 81% of BGCs lack detectable similarity to characterized MIBiG pathways
- BGC abundance and novelty show strong phylogenetic structure
- Niche-specific BGC classes: RiPP-like, ectoine, betalactone in aquatic; T3PKS in plant-associated
- GCFs show evolutionary stratification from ancient to recent acquisitions

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
| 09 | BiG-SCAPE | GCF clustering (c0.3/0.5/0.7) | 30 min |
| 10 | GCF prevalence | Family distribution across genomes | 5 min |
| 11 | Phylogenetic stats | Pagel's lambda, logistic regression | 5 min |
| 12 | Prioritize | Score BGCs for follow-up | 1 min |
| 13 | Make figures | Publication-grade figures | 1 min |
| 14 | Ecological | Niche association (Fisher exact) | 1 min |
| 15 | Phylostratigraphy | GCF evolutionary ages | 1 min |

**Total runtime: ~25-30 hours** on a 16-core, 128 GB node.

---

## Outputs

After running the pipeline:

- figures/ - 7 publication-grade figures (PDF + PNG, 300 dpi)
- results/06_antismash/parsed/ - BGC tables (663 physical regions)
- results/09_bigscape/parsed/ - GCF prevalence across 3 cutoffs
- results/11_stats/ - Pagel's lambda, logistic regression, ecological, phylostratigraphy
- results/12_prioritization/ - 7 high-priority candidate BGCs
- tables/ - manuscript tables

---

## Terminology (Important)

- **Putatively novel BGC** - no detectable KnownClusterBlast match to MIBiG v3.1. Does NOT prove novel chemistry.
- **Related BGC** - detectable MIBiG match at less than 70% gene similarity.
- **Known BGC** - detectable MIBiG match at 70% or more gene similarity. Study-specific cutoff.
- **GCF** - Gene Cluster Family
- **Rare GCF** - prevalence less than 10% of genomes
- **Core GCF** - prevalence 90% or more of genomes

See docs/08_terminology.md for full definitions.

---

## Tools and Versions

| Tool | Version | Purpose |
|------|---------|---------|
| antiSMASH | 7.1.0 | BGC prediction |
| CheckM2 | 1.1.0 | Genome quality |
| GTDB-Tk | 2.7.2 | Taxonomy |
| BiG-SCAPE | 2.0.3 | GCF clustering |
| FastTree | 2.2.0 | Phylogeny |
| R | 4.4+ | Statistics |
| MIBiG | 3.1 | Reference BGCs |
| GTDB | R232 | Reference taxonomy |
| Pfam | current | Protein domains |

---

## Reproducibility

- Every stage logs to logs/ with timestamps
- Output tables include SHA256 checksums
- Full environment captured in environment.yml
- Database versions in docs/04_databases.md

Verify any stage output:

    cd results/06_antismash/parsed
    sha256sum -c checksums.sha256

Full verification protocol: docs/10_reproducibility.md.

---

## Citation

If you use this pipeline, please cite:

> Prabhu A, Madar IH, et al. Genome-scale comparative genome mining reveals extensive biosynthetic diversity, putative novelty, and phylogenetically structured secondary metabolism in the genus Devosia. Manuscript in preparation, 2026.

See CITATION.cff for machine-readable citation info.

---

## Contact

- **Corresponding author:** Dr Inamul Hasan Madar (CIODS, Yenepoya University)
- **Pipeline maintainer:** Arvind Prabhu
- **Issues:** Please open a GitHub issue

---

## License

MIT License. See LICENSE.

---

## Environment Notes

This pipeline was developed and tested on a Linux HPC cluster with the
following directory structure:

    /dgxb_home/se26plsc001/arvind/Devosia_BGC/

Scripts assume this base path via the `base_dir` variable at the top of
each script. To run on a different system, edit `base_dir` in each script
in `scripts/stages/` and `scripts/figures/` to match your local path.

A fresh clone test with symlinked data confirmed that all figure scripts
(Fig 3-8, S1-S6) regenerate successfully from a clean checkout.
Fig 1 (BioRender) and Fig 2 (iTOL) are documented as manual workflows.

**Portability roadmap** (future work):
- Replace hardcoded paths with `Sys.getenv("DEVOSIA_ROOT")`
- Add a Docker/Singularity image for containerized reproducibility

## Data Availability

- **Code**: This GitHub repository
- **Processed data**: `results/` (small tables tracked in git)
- **Raw data**: Genome accessions listed in `data/metadata/genome_manifest.tsv`;
  regenerate via `bash run_pipeline.sh`
- **Large outputs** (antiSMASH, BiG-SCAPE, databases): excluded from git via `.gitignore`;
  regenerate via the pipeline
