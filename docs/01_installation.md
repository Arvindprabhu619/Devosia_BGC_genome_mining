# Installation

## Prerequisites

- Linux HPC with SLURM or equivalent
- 16+ cores, 128+ GB RAM
- 500 GB free disk space
- Internet access

## Step 1 - Clone the repository

    git clone https://github.com/Arvindprabhu619/Devosia_BGC_genome_mining.git
    cd Devosia_BGC_genome_mining

## Step 2 - Install micromamba (if not present)

    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
    mkdir -p ~/micromamba/bin
    mv bin/micromamba ~/micromamba/bin/
    rm -rf bin/
    export PATH=~/micromamba/bin:$PATH
    export MAMBA_ROOT_PREFIX=~/micromamba

## Step 3 - Create the conda environment

    micromamba create -f environment.yml
    micromamba activate devosia_bgc

## Step 4 - Download databases

See docs/04_databases.md. Run the helper script:

    bash scripts/utils/download_databases.sh

Downloads CheckM2 database (~3 GB), Pfam-A.hmm (~1.5 GB), and GTDB-Tk R232 (~57 GB compressed, ~94 GB extracted). MIBiG v3.1 is bundled with antiSMASH.

## Step 5 - Configure paths

Edit scripts/00_paths.sh if any paths differ on your system.

## Step 6 - Verify installation

    bash scripts/utils/check_deps.sh

## Troubleshooting

See docs/09_troubleshooting.md.
