# Devosia BGC Genome Mining

Reproducible pipeline for genome-scale BGC analysis of the genus *Devosia*.

## Pipeline stages

| # | Script | Purpose |
|---|--------|---------|
| 01 | download_genomes.sh | Fetch Devosia assemblies from NCBI |
| 02 | curate_metadata.py | Remove MAGs, duplicates, fragmented |
| 03 | build_manifest.py | Build file manifest |
| 04 | run_checkm2.sh | Genome completeness/contamination |
| 04b | parse_checkm2.py | Parse CheckM2 output |
| 05 | run_gtdbtk.sh | Taxonomic validation (GTDB R232) |
| 06 | run_antismash.sh | BGC prediction |
| 07 | extract_bgcs.py | Parse antiSMASH JSON |
| 08 | classify_novelty.py | KnownClusterBlast vs MIBiG |
| 09 | run_bigscape.sh | GCF clustering |
| 10 | gcf_prevalence.py | GCF prevalence across genomes |
| 11 | phylogenetic_stats.R | Pagel's λ, logistic regression |
| 12 | prioritize_candidates.py | Score and rank BGCs |
| 13 | make_figures.R | Publication-grade figures |

## Tools

antiSMASH 7.1.0, CheckM2 1.1.0, GTDB-Tk 2.7.2, BiG-SCAPE 2.0.3, FastTree 2.2.0, R 4.5.3

## Databases

MIBiG 3.1, Pfam current, GTDB R232, CheckM2 2021
