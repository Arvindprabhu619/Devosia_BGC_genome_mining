# Devosia BGC Genome Mining

Genome-scale comparative genome mining of the genus **Devosia**: biosynthetic diversity, ecological association, and evolutionary stratification of secondary metabolism.

[![Pipeline](https://img.shields.io/badge/pipeline-reproducible-brightgreen)](./run_pipeline.sh)
[![License](https://img.shields.io/badge/license-MIT-blue)](./LICENSE)

---

## Overview

This repository contains the complete, reproducible computational pipeline for characterizing the biosynthetic gene cluster (BGC) repertoire of the bacterial genus *Devosia* from publicly available genomes.

**Key features:**
- **14-stage pipeline** covering genome retrieval, quality control, taxonomy validation, phylogeny, BGC prediction, GCF analysis, and statistics
- **One-command reproducibility**: `bash run_pipeline.sh` regenerates all data; `bash run_figures.sh` regenerates all figures
- **8 main + 6 supplementary figures** with unified color palette
- **Complete documentation** in `docs/`

- <img width="1224" height="1285" alt="image" src="https://github.com/user-attachments/assets/fcaf793b-12e3-4de0-b310-429dedb22263" />


---

## Key Results

| Metric | Value |
|--------|-------|
| Final genomes (*Devosia* + *Devosia_A*) | 124 |
| Physical BGC regions | 663 |
| Hybrid BGC regions | 27 |
| Distinct biosynthetic classes | 21 |
| Putatively novel BGCs | 538 (81.1%) |
| Related BGCs (< 70% KCB) | 81 (12.2%) |
| Known BGCs (>= 70% KCB) | 44 (6.6%) |
| GCFs at c0.3 / c0.5 / c0.7 | 217 / 145 / 113 |
| Core GCFs (>= 90% prevalence) | 1 |
| *Devosia_A*-specific GCFs | 26 |
| Shared *Devosia*/*Devosia_A* GCFs | 14 |
| Pagel's lambda (BGC abundance) | 0.933 (P = 3.2e-11) |
| Pagel's lambda (novelty proportion) | 0.981 (P = 2.5e-29) |
| Significant niche-BGC associations (FDR < 0.05) | 5 |
| Ancient GCFs | 46 |
| Singleton GCFs | 38 |
| Prioritized candidate BGCs | 7 |

**Highlights:**
- *Devosia* genomes encode 663 BGCs across 21 biosynthetic classes
- 81% of BGCs lack detectable similarity to characterized MIBiG pathways
- BGC abundance and novelty show strong phylogenetic structure
- Niche-specific BGC classes: RiPP-like, ectoine, betalactone in aquatic; T3PKS in plant-associated
- GCFs show evolutionary stratification from ancient to recent acquisitions

---

## Repository Structure

```text
Devosia_BGC_genome_mining/
├── run_pipeline.sh
├── run_stage.sh
├── run_figures.sh
├── environment.yml
├── CITATION.cff
├── LICENSE
│
├── scripts/
│   ├── stages/
│   ├── figures/
│   └── utils/
│
├── archive/
│
├── figures/
│
├── results/
│   ├── 03_checkm2/
│   ├── 04_gtdbtk/
│   ├── 05_phylogeny/
│   ├── 06_antismash/
│   │   └── parsed/
│   ├── 09_bigscape/
│   │   └── parsed/
│   ├── 11_stats/
│   └── 12_prioritization/
│
├── itol_exports/
├── data/
│   └── metadata/
├── docs/
├── tables/
└── manuscript/
```

---

## Quick Start

### Prerequisites

- Linux HPC with SLURM or equivalent
- 16+ CPU cores
- 128+ GB RAM recommended
- 500 GB free disk space
- `micromamba` or `conda`

### Installation

```bash
git clone https://github.com/Arvindprabhu619/Devosia_BGC_genome_mining.git
cd Devosia_BGC_genome_mining
micromamba create -f environment.yml
micromamba activate devosia_bgc
bash scripts/utils/download_databases.sh
bash scripts/utils/check_deps.sh
```

### Run the Full Pipeline

```bash
bash run_pipeline.sh
```

### Run a Range of Stages

```bash
bash run_pipeline.sh 06 10
```

### Run a Single Stage

```bash
bash run_stage.sh 06
```

### Regenerate Figures

```bash
bash run_figures.sh
```

---

## Pipeline Stages

| Stage | Name | Purpose | Approx. Runtime |
|------:|------|---------|-----------------|
| 01 | Download genomes | Fetch *Devosia* assemblies from NCBI | 15–30 min |
| 02 | Curate metadata | Remove MAGs, duplicates, fragmented genomes | 2 min |
| 03 | Build manifest | Map accessions to genome files | 2 min |
| 04 | CheckM2 | Genome completeness/contamination | 20–60 min |
| 05 | GTDB-Tk classify | Taxonomic validation using GTDB R232 | 20–60 min |
| 05b | Build phylogeny | bac120 alignment + FastTree | 30–60 min |
| 06 | antiSMASH | BGC prediction | 20–24 hr |
| 07 | Extract BGCs | Parse antiSMASH JSON | 5 min |
| 08 | Classify novelty | KnownClusterBlast against MIBiG | 5 min |
| 09 | BiG-SCAPE | GCF clustering | 30 min |
| 10 | GCF prevalence | Family distribution | 5 min |
| 11 | Phylogenetic statistics | Pagel's lambda, logistic regression | 5 min |
| 12 | Prioritize | Score BGCs for follow-up | 1 min |
| 13 | Ecological | Niche association | 1 min |
| 14 | Phylostratigraphy | GCF evolutionary ages | 1 min |

**Note:** Figures are generated separately using `run_figures.sh`.

---

## Figure List

### Main Figures

| Figure | Script / Workflow | Description |
|--------|-------------------|-------------|
| Fig 1 | BioRender (manual) | Study workflow - 516 to 124 genome curation |
| Fig 2 | iTOL (manual) | Circular phylogeny with 20 BGC class rings |
| Fig 3 | `03_biosynthetic_diversity.R` | Biosynthetic composition, heatmap, and novelty |
| Fig 4 | `04_gcf_diversity.R` | GCF prevalence, ranked abundance, and overlap |
| Fig 5 | `05_phylo_structure.R` | Pagel's lambda and logistic regression |
| Fig 6 | `06_ecological_evolutionary.R` | Niche associations and GCF ages |
| Fig 7 | `07_prioritized_candidates.R` | Candidate BGC classes and prioritization scores |
| Fig 8 | `08_bgc_neighborhoods.R` | Representative BGC gene neighborhoods |

### Supplementary Figures

| Figure | Script / Workflow | Description |
|--------|-------------------|-------------|
| S1 | iTOL (manual) | Multitrack circular tree |
| S2 | `S2_checkm2_qc.R` | CheckM2 quality assessment |
| S3 | `S3_candidate_details.R` | Candidate BGC details |
| S4 | `S4_novelty_by_genus.R` | Novelty by genus |
| S5 | `S5_rarefaction.R` | GCF rarefaction |
| S6 | `S6_bgc_lengths.R` | BGC length distribution |

Fig 1 and Fig 2 are manual workflows.

Fig 2 inputs are in `results/itOL_exports/`

```text
results/itOL_exports/
```

---

## Terminology

| Term | Definition |
|------|------------|
| Putatively novel BGC | No KnownClusterBlast match to MIBiG v3.1 |
| Related BGC | MIBiG match < 70% gene similarity |
| Known BGC | MIBiG match >= 70% gene similarity |
| GCF | Gene Cluster Family |
| Rare GCF | < 10% of genomes |
| Core GCF | >= 90% of genomes |

See `docs/08_terminology.md` for full definitions.

---

## Tools and Versions

| Tool | Version | Purpose |
|------|---------|---------|
| antiSMASH | 7.1.0 | BGC prediction |
| CheckM2 | 1.1.0 | Genome quality |
| GTDB-Tk | 2.7.2 | Taxonomy |
| BiG-SCAPE | 2.0.3 | GCF clustering |
| FastTree | 2.2.0 | Phylogeny |
| R | 4.1.2 | Statistics and figures |
| ggplot2 | 3.5.1 | Visualization |
| patchwork | 1.2.0 | Multi-panel figure assembly |
| MIBiG | 3.1 | Reference BGCs |
| GTDB | R232 | Reference taxonomy |
| Pfam | current | Protein domains |

---

## Reproducibility

Every pipeline stage logs execution information to:

```text
logs/
```

SHA256 checksums are provided for selected processed outputs.

Environment configuration is defined in:

```text
environment.yml
```

Database versions and acquisition information are documented in:

```text
docs/04_databases.md
```

### Verify Output

```bash
cd results/06_antismash/parsed
sha256sum -c checksums.sha256
```

Figures 3–8 and supplementary figures can be regenerated using:

```bash
bash run_figures.sh
```

Figures 1–2 are manual workflows.

---

## Expected Output

The main processed results are organized by analytical stage:

```text
results/
├── 03_checkm2/
├── 04_gtdbtk/
├── 05_phylogeny/
├── 06_antismash/
│   └── parsed/
├── 09_bigscape/
│   └── parsed/
├── 11_stats/
├── 12_prioritization/
├── 13_ecology/
└── 14_phylostratigraphy/
```

Large intermediate files and computationally intensive outputs are excluded from version control where appropriate and can be regenerated through the pipeline.

---

## Study Summary

The final curated dataset contains **124 genomes** representing *Devosia* and *Devosia_A*. Genome-scale mining identified **663 physical BGC regions** spanning **21 biosynthetic classes**.

Of these, **538 BGCs (81.1%)** were classified as putatively novel based on the absence of detectable KnownClusterBlast similarity to MIBiG v3.1 under the defined thresholding framework.

A further **81 BGCs (12.2%)** showed relatedness below the 70% gene-similarity threshold, while **44 BGCs (6.6%)** met the defined criterion for known BGCs.

GCF analysis identified:

- 217 GCFs at c0.3
- 145 GCFs at c0.5
- 113 GCFs at c0.7
- 1 core GCF
- 26 *Devosia_A*-specific GCFs
- 14 GCFs shared between *Devosia* and *Devosia_A*

Phylogenetic analyses identified significant phylogenetic structure in BGC abundance and BGC novelty.

Ecological analyses identified significant associations between selected biosynthetic classes and ecological niches.

Phylostratigraphic analysis categorized GCFs according to their inferred evolutionary distribution.

Seven candidate BGCs were prioritized for potential downstream investigation using the defined prioritization framework.

---

## Citation

Prabhu A, Madar IH, et al. **Genome-scale comparative genome mining reveals extensive biosynthetic diversity, putative novelty, and phylogenetically structured secondary metabolism in the genus *Devosia*.** Manuscript in preparation, 2026.

See `CITATION.cff` for machine-readable citation metadata.

---

## Contact

**Corresponding author:**  
Dr Inamul Hasan Madar  
Centre for Integrative Omics and Data Sciences (CIODS)  
Yenepoya University

**Pipeline maintainer:**  
Arvind Prabhu

For questions, reproducibility issues, or proposed improvements, please open a GitHub issue.

---

## License

MIT License.

See `LICENSE` for the complete license text.

---

## Environment Notes

This pipeline was tested on a Linux HPC cluster with the following project directory:

```text
/dgxb_home/se26plsc001/arvind/Devosia_BGC/
```

Scripts currently assume this base path through a `base_dir` variable.

To run the pipeline elsewhere, update the relevant `base_dir` settings in:

```text
scripts/stages/
scripts/figures/
```

### HPC Requirements

Recommended resources:

- >=16 CPU cores
- >=128 GB RAM
- >=500 GB available storage
- Linux operating system
- Conda or micromamba

### Fresh-Clone Validation

A fresh-clone test confirmed that the automated figure-generation workflow can regenerate:

- Fig 3–8
- S1–S6

The following remain manual workflows:

- Fig 1 — BioRender
- Fig 2 — iTOL

---

## Portability Roadmap

Future improvements include:

1. Replace hard-coded paths with `DEVOSIA_ROOT`.

For example:

```bash
export DEVOSIA_ROOT=/path/to/Devosia_BGC
```

and in R scripts:

```r
Sys.getenv("DEVOSIA_ROOT")
```

2. Add Docker or Apptainer/Singularity support.

3. Further automate iTOL annotation-file generation.

4. Add automated validation of expected output files after each pipeline stage.

5. Add continuous-integration testing for core scripts.

---

## Data Availability

### Code

The complete analysis code is available in this GitHub repository.

### Processed Data

Processed tables and summary outputs are available under:

```text
results/
tables/
```

### Raw Genome Data

Raw genome files are not stored directly in the repository.

Genome accession information is maintained in:

```text
data/metadata/genome_manifest.tsv
```

The genome dataset can be regenerated using:

```bash
bash run_pipeline.sh
```

### Large Outputs

Large intermediate and computational outputs are excluded from Git version control using `.gitignore` where appropriate.

These files can be regenerated through the documented pipeline.

---

## Reproducibility Checklist

Before considering a complete reproduction of the analysis, verify:

```text
[ ] environment.yml successfully creates the analysis environment
[ ] Required databases are downloaded
[ ] Genome manifest is generated
[ ] CheckM2 results are generated
[ ] GTDB-Tk classification is completed
[ ] Phylogenetic tree is generated
[ ] antiSMASH results are available
[ ] BGC parsing is completed
[ ] MIBiG novelty classification is completed
[ ] BiG-SCAPE GCF analysis is completed
[ ] GCF prevalence is calculated
[ ] Phylogenetic statistics are completed
[ ] Ecological analysis is completed
[ ] Phylostratigraphy is completed
[ ] Candidate prioritization is completed
[ ] Figures 3–8 are regenerated
[ ] Supplementary figures are regenerated
[ ] iTOL inputs are available for manual figures
```

---

## One-Command Reproduction

After installation and database preparation:

```bash
cd ~/arvind/Devosia_BGC/
bash run_pipeline.sh
```

To regenerate computational figures:

```bash
bash run_figures.sh
```

To inspect repository status:

```bash
git status
```

---

## Project Status

**Status:** Manuscript preparation / computational analysis complete

The repository contains the computational workflow, processed analysis outputs, figure-generation scripts, documentation, and manuscript-supporting tables used for the comparative genome-mining study of *Devosia*.

---

## Acknowledgements

We acknowledge the developers and maintainers of NCBI, antiSMASH, CheckM2, GTDB-Tk, BiG-SCAPE, FastTree, MIBiG, Pfam, R, and associated open-source software and databases used in this study.

We also acknowledge the public genome repositories and sequencing projects that made the comparative analysis possible.

---

## Repository

GitHub:

https://github.com/Arvindprabhu619/Devosia_BGC_genome_mining
