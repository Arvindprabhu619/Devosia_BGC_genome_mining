# Expected Outputs

After the pipeline completes, you will find outputs in these locations.

## results/

### results/04_checkm2/
- quality_report.tsv - per-genome completeness/contamination
- quality_passed.txt - accessions passing QC
- quality_excluded.txt - accessions failing QC

### results/05_gtdbtk/
- final_genome_list.txt - Devosia/Devosia_A accessions
- taxonomy_assignments.tsv - full taxonomy per genome
- excluded_taxa.txt - non-target genera

### results/05b_phylogeny/
- devosia_pruned.tree - final phylogenetic tree (Newick)
- devosia_alignment.fasta - bac120 MSA

### results/06_antismash/
- raw_output/<acc>/<acc>.gbk - GenBank format BGCs
- raw_output/<acc>/<acc>.json - JSON format BGCs
- parsed/bgc_regions.tsv - physical BGC regions
- parsed/bgc_class_assignments.tsv - class-expanded BGCs
- parsed/bgc_novelty.tsv - per-BGC novelty classification

### results/09_bigscape/
- parsed/gcf_prevalence_c0.3.tsv
- parsed/gcf_prevalence_c0.5.tsv
- parsed/gcf_prevalence_c0.7.tsv
- parsed/gcf_overlap_devosia_A.tsv

### results/11_stats/
- pagels_lambda.tsv - phylogenetic signal
- logistic_regression.tsv - class associations
- ecological_association_classes.tsv
- ecological_association_gcfs.tsv
- gcf_ages.tsv - phylostratigraphy

### results/12_prioritization/
- scored_bgcs.tsv - all scored BGCs
- candidate_pool.tsv - BGCs above threshold
- prioritized_candidates.tsv - top 15

## figures/

- Fig2_genome_landscape.pdf / .png
- Fig3_bgc_classes.pdf / .png
- Fig4_novelty_landscape.pdf / .png
- Fig5_phylogenetic_structure.pdf / .png
- Fig6_gcf_analysis.pdf / .png
- Fig7_ecological_phylostrat.pdf / .png
- Fig8_prioritized_candidates.pdf / .png

All figures at 300 dpi.

## logs/

Per-stage logs with timestamps. Format: <stage>_<date>_<time>.log
