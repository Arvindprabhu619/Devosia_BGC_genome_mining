# Databases

## Required databases

| Database | Version | Size | Purpose |
|----------|---------|------|---------|
| CheckM2 DB | 2021 release | ~3 GB | Genome quality |
| Pfam-A.hmm | current | ~1.5 GB | Protein domains |
| GTDB-Tk R232 | R232 | 57 GB tar | Taxonomy + markers |
| MIBiG | 3.1 | bundled | Reference BGCs |

## Download

    bash scripts/utils/download_databases.sh

## Manual downloads

### CheckM2 DB

    cd data/databases
    mkdir -p checkm2_db && cd checkm2_db
    wget https://zenodo.org/records/5571251/files/checkm2_database.tar.gz
    tar -xzf checkm2_database.tar.gz

### Pfam-A.hmm

    cd data/databases
    wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm

### GTDB-Tk R232

    cd data/databases
    mkdir -p gtdbtk_r232 && cd gtdbtk_r232
    wget https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
    tar -xzf gtdbtk_r232_data.tar.gz

### MIBiG 3.1

Bundled with antiSMASH 7.1.0. No separate download needed.

## Verification

    ls -lh data/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd
    ls -lh data/databases/Pfam-A.hmm*
    ls data/databases/gtdbtk_r232/release232/
