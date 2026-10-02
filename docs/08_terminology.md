# Terminology

## BGC-related terms

**Biosynthetic Gene Cluster (BGC)** - Genomic region containing functionally associated genes involved in biosynthesis, modification, regulation, transport, and protection of a specialized metabolite.

**Physical BGC region** - One antiSMASH-predicted region. Counted once regardless of number of product classes.

**Hybrid BGC** - A physical region carrying >=2 biosynthetic product classifications.

**Class assignment** - A single biosynthetic class from antiSMASH. Hybrid BGCs contribute multiple class assignments.

## Novelty classification (IMPORTANT)

**Putatively novel BGC** - A BGC with NO detectable KnownClusterBlast match to MIBiG v3.1. Does NOT prove novel chemistry. Indicates absence from the current MIBiG reference collection under the analysis conditions.

**Related BGC** - A BGC with a detectable KnownClusterBlast match at < 70% gene similarity. Shows sequence-level relatedness but divergence from the characterized reference.

**Known BGC** - A BGC with a detectable KnownClusterBlast match at >= 70% gene similarity. Study-specific operational cutoff. Does NOT imply confirmed production of the reference compound.

**The 70% threshold** - A study-specific operational boundary, not a universal definition of BGC identity or functional equivalence. antiSMASH KnownClusterBlast reports the proportion of reference genes with detectable sequence-similar matches in the query.

## Gene cluster family terms

**GCF (Gene Cluster Family)** - A cluster of BGCs with similar domain content and organization, defined by BiG-SCAPE.

**Rare GCF** - Prevalence < 10% of analyzed genomes.

**Intermediate GCF** - Prevalence 10-90%.

**Core GCF** - Prevalence >= 90%.

**Singleton GCF** - A GCF with only one member BGC (retained via --include-singletons).

## Phylogenetic terms

**Pagel lambda** - A measure of phylogenetic signal. lambda = 0 indicates no phylogenetic dependence; lambda = 1 indicates trait covariance follows Brownian expectations.

**Phylogenetic logistic regression** - Regression of binary traits (presence/absence) accounting for phylogenetic covariance.

**MRCA** - Most Recent Common Ancestor.

**Phylostratigraphy** - Assigning each GCF to an evolutionary age based on MRCA depth.

## Ecological terms

**Niche categories** - soil, plant, aquatic, contaminated, host-associated, industrial. Assigned via keyword patterns in BioSample metadata.

**Fisher exact test** - Statistical test for 2x2 contingency tables, used to test niche x BGC associations.
