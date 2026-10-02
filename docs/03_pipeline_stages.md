# Pipeline Stages

Detailed documentation for each stage.

## Stage 01 - Download Genomes

Script: scripts/01_download_genomes.sh
Tool: NCBI Datasets CLI 18.38.0
Runtime: 15-30 min

Queries NCBI for all *Devosia* assemblies. Outputs data/metadata/devosia_refseq.jsonl, devosia_genbank.jsonl, and data/raw_genomes/ncbi_dataset/.

## Stage 02 - Curate Metadata

Script: scripts/02_curate_metadata.py. Runtime: 2 min.

Removes MAGs, uncultured/Candidatus, and fragmented assemblies (>200 scaffolds). Consolidates GenBank/RefSeq duplicates.

## Stage 03 - Build Manifest

Script: scripts/03_build_manifest.sh. Runtime: 2 min.

Maps each accession to its FASTA/GFF3/FAA files.

## Stage 04 - CheckM2

Script: scripts/04_run_checkm2.sh. Tool: CheckM2 1.1.0. Runtime: 20-60 min.

Thresholds: standard = completeness >=95% and contamination <=5.0%; borderline = completeness >=95% and contamination 5.0-5.5%; excluded otherwise.

## Stage 05 - GTDB-Tk

Script: scripts/05_run_gtdbtk.sh. Tool: GTDB-Tk 2.7.2. Runtime: 20-60 min.

Assigns each genome to GTDB R232 taxonomy. Symlinks prefixed with query_ to avoid ID collision. Retains only Devosia and Devosia_A.

## Stage 05b - Phylogeny

Script: scripts/05b_build_phylogeny.sh. Tools: GTDB-Tk + FastTree 2.2.0. Runtime: 30-60 min.

Two-step: gtdbtk identify (bac120 markers) then gtdbtk align; FastTree with WAG+Gamma. Output: devosia_pruned.tree.

## Stage 06 - antiSMASH

Script: scripts/06_run_antismash.sh. Tool: antiSMASH 7.1.0. Runtime: 20-24 hr.

Predicts BGCs. 4 genomes processed in parallel. Outputs gbk + json files per genome.

## Stage 07 - Extract BGCs

Script: scripts/07_extract_bgcs.py. Runtime: 5 min.

Parses antiSMASH JSON into bgc_regions.tsv, bgc_class_assignments.tsv, hybrid_regions.tsv.

## Stage 08 - Classify Novelty

Script: scripts/08_classify_novelty.py. Runtime: 5 min.

Compares each BGC to MIBiG via KnownClusterBlast. Categories: putatively novel (no match), related (<70%), known (>=70%).

## Stage 09 - BiG-SCAPE

Script: scripts/09_run_bigscape.sh. Tool: BiG-SCAPE 2.0.3. Runtime: 2-4 hr.

Clusters BGCs into GCFs at cutoffs 0.3, 0.5, 0.7.

## Stage 10 - GCF Prevalence

Script: scripts/10_gcf_prevalence.py. Runtime: 15 min.

Categories: core (>=90%), intermediate (10-90%), rare (<10%).

## Stage 11 - Phylogenetic Stats

Script: scripts/11_phylogenetic_stats.R. Runtime: 1-2 hr.

Pagel lambda for BGC abundance and novelty proportion; phylogenetic logistic regression per BGC class.

## Stage 12 - Prioritize Candidates

Script: scripts/12_prioritize_candidates.py. Runtime: 15 min.

Multi-criterion scoring. Selects top 15.

## Stage 13 - Figures

Script: scripts/13_make_figures.R. Runtime: 30 min.

Publication-grade figures (PDF + PNG, 300 dpi).

## Stage 14 - Ecological Association

Script: scripts/14_ecological_association.py. Runtime: 30 min.

Fisher exact test for niche x class / niche x GCF.

## Stage 15 - Phylostratigraphy

Script: scripts/15_phylostratigraphy.R. Runtime: 1 hr.

Assigns each GCF an evolutionary age via MRCA depth.
