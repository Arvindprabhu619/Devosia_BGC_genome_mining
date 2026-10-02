#!/usr/bin/env bash
# ============================================================
# scripts/utils/check_deps.sh
# ============================================================
# Purpose: Verify all pipeline dependencies are installed.
# Usage:   bash scripts/utils/check_deps.sh
# ============================================================
set -uo pipefail
cd "$(dirname "$0")/../.."
source scripts/00_paths.sh

echo "=== Devosia BGC pipeline dependency check ==="
echo ""

echo "--- System ---"
echo "Hostname:  $(hostname)"
echo "CPUs:      $(nproc)"
echo "Free disk: $(df -h ~ | tail -1 | awk "{print \$4}")"
echo ""

echo "--- Tool wrappers ---"
MISSING=0
for tool in ANTI_SMASH CHECKM2 GTDBTK BIGSCAPE FASTTREE DATASETS; do
    wrapper_var=$(eval echo "\$$tool")
    if [[ -x "$wrapper_var" ]]; then
        printf "%-15s OK\n" "$tool"
    else
        printf "%-15s MISSING (%s)\n" "$tool" "$wrapper_var"
        MISSING=$((MISSING+1))
    fi
done
echo ""

echo "--- R wrapper ---"
if [[ -x "$RSCRIPT" ]]; then
    echo "Rscript: OK"
else
    echo "Rscript: MISSING"
    MISSING=$((MISSING+1))
fi
echo ""

echo "--- Python ---"
if [[ -x "$PYTHON_BIN" ]]; then
    echo "Python: $PYTHON_BIN"
    "$PYTHON_BIN" --version
else
    echo "Python: MISSING ($PYTHON_BIN)"
    MISSING=$((MISSING+1))
fi
echo ""

echo "--- Databases ---"
for db in CHECKM2_DB PFAM_HMM GTDBTK_DATA; do
    db_path=$(eval echo "\$$db")
    if [[ -e "$db_path" ]]; then
        printf "%-15s OK (%s)\n" "$db" "$(du -sh "$db_path" 2>/dev/null | cut -f1)"
    else
        printf "%-15s MISSING (%s)\n" "$db" "$db_path"
        MISSING=$((MISSING+1))
    fi
done
echo ""

echo "--- Directories ---"
for d in DATA_RAW DATA_META DATA_DB RES_ANTISMASH DIR_FIGS; do
    d_path=$(eval echo "\$$d")
    if [[ -d "$d_path" ]]; then
        printf "%-15s OK\n" "$d"
    else
        printf "%-15s MISSING\n" "$d"
    fi
done
echo ""

if [[ $MISSING -eq 0 ]]; then
    echo "All dependencies OK"
    exit 0
else
    echo "$MISSING dependencies MISSING - see docs/09_troubleshooting.md"
    exit 1
fi
