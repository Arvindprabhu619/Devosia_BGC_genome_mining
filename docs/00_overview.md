# Overview

## Scientific question

The genus *Devosia* (family Devosiaceae, order Hyphomicrobiales) comprises environmentally diverse alphaproteobacteria found in soils, rhizospheres, freshwater, marine habitats, and contaminated environments. While individual *Devosia* species have been studied for metabolic versatility, the genus-wide biosynthetic potential for secondary metabolites has not been systematically characterized.

## Hypothesis

The broad ecological and genomic diversity of *Devosia* is reflected in a heterogeneous biosynthetic repertoire containing: conserved BGC families across many genomes, genome-restricted accessory BGC families, substantial sequence space lacking similarity to characterized pathways, BGCs associated with specific ecological niches, and BGCs stratified by evolutionary age.

## Approach

1. Curate publicly available *Devosia* genomes with quality control
2. Validate taxonomy using GTDB-Tk R232
3. Reconstruct a genome-scale phylogeny
4. Predict BGCs with antiSMASH 7.1.0
5. Assess novelty against MIBiG v3.1 via KnownClusterBlast
6. Cluster BGCs into gene cluster families (BiG-SCAPE)
7. Test phylogenetic structuring (Pagel lambda, logistic regression)
8. Associate BGCs with ecological niches (Fisher exact)
9. Date GCFs via phylostratigraphy
10. Prioritize candidate BGCs for experimental follow-up

## Novel contributions

- **Ecological structure**: BGC classes and families are non-randomly distributed across isolation niches
- **Evolutionary stratification**: GCFs span a range of evolutionary ages from ancient to recent
- **Reproducible workflow**: full pipeline in one repository, one command
