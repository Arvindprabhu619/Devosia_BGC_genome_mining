
# Genome-wide mining of biosynthetic gene clusters in *Devosia*

## Overview

This repository contains the computational workflows, curated metadata,
analysis results, phylogenetic data, and publication outputs generated for
the genome-wide analysis of biosynthetic gene clusters (BGCs) in the genus
*Devosia*.

The analysis integrates genome curation, genome-quality assessment,
taxonomic validation, genome statistics, phylogenetic reconstruction,
antiSMASH BGC prediction, BiG-SCAPE gene-cluster-family (GCF) analysis,
phylogenetic statistics, and candidate BGC prioritization.

The repository provides a transparent and reproducible record of the
computational analyses supporting the associated manuscript.

---

## Dataset

The genome dataset was assembled from NCBI assembly records and underwent
successive curation and quality-control steps.

### Genome-selection workflow

1. NCBI assembly records were retrieved and consolidated.
2. Paired GenBank (GCA) and RefSeq (GCF) assemblies representing the same
   physical assembly were consolidated.
3. Cultured assemblies were retained for downstream analysis.
4. Genome quality was evaluated using CheckM2.
5. Taxonomic identity was assessed using GTDB-Tk.
6. Non-*Devosia* genomes were excluded.
7. The final dataset contained 110 taxonomically validated genomes.

### Final dataset

- Total genomes: **110**
- *Devosia*: **88**
- *Devosia_A*: **22**
- Final phylogeny: **110 tips**

---

## Computational workflow

The repository is organized according to the major stages of the analysis:

```text
01_data_curation
        |
        v
02_quality_control
        |
        v
03_taxonomy
        |
        v
04_genome_statistics
        |
        v
05_phylogeny
        |
        v
06_antismash
        |
        v
07_bigscape
        |
        v
08_gcf_analysis
        |
        v
09_phylogenetic_statistics
        |
        v
10_candidate_prioritization
        |
        v
11_visualization
