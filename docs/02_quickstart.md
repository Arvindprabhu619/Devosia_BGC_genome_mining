# Quickstart

## One-command run

    bash run_pipeline.sh

Runs all 15 stages sequentially. Estimated runtime: 30-40 hours on a 16-core node.

## Stage-by-stage

    bash run_pipeline.sh 01 05      # stages 01-05
    bash run_pipeline.sh 06 06      # stage 06 only
    bash run_stage.sh 06            # single stage (debugging)
    bash run_stage.sh               # list all stages

## What you get

After the pipeline completes:

- figures/ : 7 publication-grade PDF + PNG figures
- tables/ : manuscript tables
- results/12_prioritization/prioritized_candidates.tsv : top 15 BGCs
- results/11_stats/ : phylogenetic + ecological statistics
- logs/ : every stage log

## Verify

    ls results/06_antismash/raw_output | wc -l

## Resuming after failure

All scripts are idempotent. Fix the issue, then:

    bash run_pipeline.sh 06 15
