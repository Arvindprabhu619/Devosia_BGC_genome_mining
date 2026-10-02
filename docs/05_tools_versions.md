# Tools and Versions

## Core tools

| Tool | Version | Purpose |
|------|---------|---------|
| antiSMASH | 7.1.0 | BGC prediction |
| CheckM2 | 1.1.0 | Genome completeness/contamination |
| GTDB-Tk | 2.7.2 | Taxonomic classification |
| BiG-SCAPE | 2.0.3 | GCF clustering |
| FastTree | 2.2.0 | Phylogeny |
| Prodigal | 2.6.3 | Gene prediction |
| HMMER | 3.4 | Profile HMM searches |
| DIAMOND | 2.x | Protein alignment |
| NCBI Datasets CLI | 18.38.0 | Genome download |

## R packages

| Package | Version | Purpose |
|---------|---------|---------|
| ape | 5.x | Phylogenetics |
| phytools | 2.x | Comparative methods |
| phylolm | 2.x | Phylogenetic regression |
| ggtree | 3.x | Tree visualization |
| ggplot2 | 3.x | Plotting |
| dplyr, tidyr | latest | Data manipulation |
| patchwork, cowplot | latest | Multi-panel figures |
| viridis, RColorBrewer | latest | Color palettes |

## Python packages

- pandas, numpy, scipy, statsmodels, biopython

## Reference databases

| Database | Version |
|----------|---------|
| MIBiG | 3.1 |
| GTDB | R232 |
| Pfam | current |

## Capture versions

    bash scripts/utils/check_deps.sh > env/software_versions.txt
