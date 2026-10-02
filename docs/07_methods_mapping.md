# Manuscript Methods to Code Mapping

Each Methods subsection in the manuscript corresponds to specific scripts and parameters.

## 2.1 Genome retrieval and metadata curation

Scripts: scripts/01_download_genomes.sh, scripts/02_curate_metadata.py

- NCBI query: taxon "Devosia", RefSeq + GenBank
- Exclude MAGs: keyword regex on organism name + BioSample metadata
- Exclude uncultured/Candidatus: keyword regex
- Exclude fragmented: > 200 scaffolds
- Consolidate GenBank/RefSeq: prefer RefSeq (GCF_)

## 2.2 Genome quality assessment and taxonomic validation

Scripts: scripts/04_run_checkm2.sh, scripts/05_run_gtdbtk.sh

- CheckM2 with default lineage-independent workflow
- Standard: completeness >= 95%, contamination <= 5.0%
- Borderline: completeness >= 95%, contamination 5.0-5.5%
- Excluded: completeness < 95% OR contamination > 5.5%
- GTDB-Tk classify_wf, GTDB R232
- Retain only Devosia and Devosia_A

## 2.4 Phylogenomic reconstruction

Script: scripts/05b_build_phylogeny.sh

- GTDB-Tk identify: bac120 marker set
- GTDB-Tk align: concatenated amino acid alignment
- FastTree: WAG + Gamma

## 2.5 Biosynthetic gene cluster prediction

Script: scripts/06_run_antismash.sh

- antiSMASH 7.1.0 with --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --rre
- Prodigal gene prediction
- Min length 1000 bp

## 2.6 BGC novelty classification

Script: scripts/08_classify_novelty.py

- KnownClusterBlast vs MIBiG 3.1
- Putatively novel: no ranked match
- Related: match < 70% gene similarity
- Known: match >= 70% gene similarity

## 2.9 Phylogenetic statistical analysis

Script: scripts/11_phylogenetic_stats.R

- Pagel lambda via phytools::phylosig
- Phylogenetic logistic regression via phylolm
- FDR correction (Benjamini-Hochberg)

## 2.10 Gene cluster family analysis

Scripts: scripts/09_run_bigscape.sh, scripts/10_gcf_prevalence.py

- BiG-SCAPE 2.0.3 with Pfam-A.hmm and MIBiG 3.1
- Cutoffs: 0.3, 0.5, 0.7
- --include-singletons
- Prevalence: core >= 90%, intermediate 10-90%, rare < 10%

## 2.12 Prioritization

Script: scripts/12_prioritize_candidates.py

- Score: novelty + rarity + class + hybrid + domains + candidates
- Threshold: score >= 6
- Top 15 selected
