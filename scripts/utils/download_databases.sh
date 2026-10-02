#!/usr/bin/env bash
# ============================================================
# scripts/utils/download_databases.sh
# ============================================================
# Purpose: Download all required databases.
# Usage:   bash scripts/utils/download_databases.sh
# ============================================================
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/00_paths.sh

LOG="${DIR_LOGS}/database_download_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== Database download: $(date) ==="
echo ""

# ---- CheckM2 DB ----
if [[ -f "$CHECKM2_DB" ]]; then
    echo "[CheckM2 DB] Already exists at $CHECKM2_DB"
else
    echo "[CheckM2 DB] Downloading..."
    mkdir -p "$(dirname "$CHECKM2_DB")/.."
    cd "$(dirname "$CHECKM2_DB")/.."
    wget -c https://zenodo.org/records/5571251/files/checkm2_database.tar.gz
    tar -xzf checkm2_database.tar.gz
    rm -f checkm2_database.tar.gz
    cd - > /dev/null
    echo "[CheckM2 DB] Done"
fi
echo ""

# ---- Pfam ----
if [[ -f "$PFAM_HMM" && -f "${PFAM_HMM}.h3f" ]]; then
    echo "[Pfam] Already exists at $PFAM_HMM"
else
    echo "[Pfam] Downloading (~1.5 GB)..."
    cd "$DATA_DB"
    wget -c http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip -f Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
    rm -f Pfam-A.hmm.gz
    cd - > /dev/null
    echo "[Pfam] Done"
fi
echo ""

# ---- GTDB-Tk R232 ----
if [[ -d "$GTDBTK_DATA" && -d "${GTDBTK_DATA}/metadata" ]]; then
    echo "[GTDB-Tk] Already extracted at $GTDBTK_DATA"
else
    echo "[GTDB-Tk] Downloading (~57 GB compressed, ~94 GB extracted)..."
    echo "[GTDB-Tk] This takes 2-4 hours."
    mkdir -p "$(dirname "$GTDBTK_DATA")"
    cd "$(dirname "$GTDBTK_DATA")"
    if [[ ! -f gtdbtk_r232_data.tar.gz ]]; then
        wget -c https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
    fi
    echo "[GTDB-Tk] Extracting..."
    tar -xzf gtdbtk_r232_data.tar.gz
    rm -f gtdbtk_r232_data.tar.gz
    cd - > /dev/null
    echo "[GTDB-Tk] Done"
fi
echo ""

# ---- MIBiG ----
echo "[MIBiG] Bundled with antiSMASH 7.1.0 - no separate download"
echo ""

# ---- Checksums ----
echo "=== Database checksums ==="
sha256sum "$CHECKM2_DB" > "${DIR_ENV}/database_versions.txt" 2>/dev/null || true
echo "Written to: ${DIR_ENV}/database_versions.txt"

echo "=== Database download complete: $(date) ==="
