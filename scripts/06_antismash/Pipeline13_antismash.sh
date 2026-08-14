#!/usr/bin/env bash

# ============================================================
# Pipeline 13 — antiSMASH BGC prediction
#
# Purpose:
#   Run antiSMASH independently for each of the final curated
#   Devosia genomes.
#
# Input:
#   final_validation/final_genomes/*.fna
#
# Output:
#   antismash_output/<accession>/
#
# Portability:
#   Repository root is inferred from this script's location.
#   No hard-coded HOME or project-specific paths.
# ============================================================

set -euo pipefail


# ------------------------------------------------------------
# 1. Determine repository root
# ------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"


# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------

FINAL_DIR="$PROJECT_ROOT/final_validation"

GENOME_DIR="$FINAL_DIR/final_genomes"

ANTISMASH_DIR="$PROJECT_ROOT/antismash_output"

LOG_DIR="$PROJECT_ROOT/logs"

LIST_FILE="$PROJECT_ROOT/antismash_110_list.tsv"

RUNNER="$PROJECT_ROOT/run_one_antismash.sh"


# ------------------------------------------------------------
# 3. Input validation
# ------------------------------------------------------------

if [ ! -d "$GENOME_DIR" ]; then
    echo "ERROR: Final genome directory not found:"
    echo "  $GENOME_DIR"
    exit 1
fi

if [ ! -f "$LIST_FILE" ]; then
    echo "ERROR: antiSMASH genome list not found:"
    echo "  $LIST_FILE"
    echo
    echo "Expected format:"
    echo "  genome_name<TAB>accession"
    exit 1
fi


# ------------------------------------------------------------
# 4. Create output directories
# ------------------------------------------------------------

mkdir -p "$ANTISMASH_DIR"
mkdir -p "$LOG_DIR"


# ------------------------------------------------------------
# 5. Resource report
# ------------------------------------------------------------

echo "============================================================"
echo "Pipeline 13 — antiSMASH BGC prediction"
echo "============================================================"

echo
echo "Project root:"
echo "  $PROJECT_ROOT"

echo
echo "Genome directory:"
echo "  $GENOME_DIR"

echo
echo "antiSMASH output:"
echo "  $ANTISMASH_DIR"

echo
echo "Log directory:"
echo "  $LOG_DIR"

echo
echo "Genome list:"
echo "  $LIST_FILE"

echo
echo "System memory:"
free -h

echo
echo "Filesystem:"
df -h "$PROJECT_ROOT"


# ------------------------------------------------------------
# 6. Create portable per-genome runner
# ------------------------------------------------------------

cat > "$RUNNER" <<'RUNNER_EOF'
#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

NAME="$1"
ACC="$2"

GENOME_DIR="$PROJECT_ROOT/final_validation/final_genomes"
ANTISMASH_DIR="$PROJECT_ROOT/antismash_output"
LOG_DIR="$PROJECT_ROOT/logs"

FASTA="$GENOME_DIR/${NAME}.fna"
OUTDIR="$ANTISMASH_DIR/${ACC}"
LOG="$LOG_DIR/${ACC}.log"

mkdir -p "$ANTISMASH_DIR"
mkdir -p "$LOG_DIR"


# ------------------------------------------------------------
# Input check
# ------------------------------------------------------------

if [ ! -f "$FASTA" ]; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: FASTA not found: $FASTA" \
        | tee -a "$LOG_DIR/master.log"
    exit 1
fi


# ------------------------------------------------------------
# Skip completed genomes
# ------------------------------------------------------------

if [ -f "${OUTDIR}/${ACC}.json" ]; then

    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ${ACC} already done, skipping" \
        | tee -a "$LOG_DIR/master.log"

    exit 0
fi


# ------------------------------------------------------------
# Clean incomplete output
# ------------------------------------------------------------

rm -rf "$OUTDIR"

START=$(date +%s)

echo "[$(date +'%Y-%m-%d %H:%M:%S')] START ${ACC}" \
    | tee -a "$LOG_DIR/master.log"


# ------------------------------------------------------------
# Run antiSMASH
# ------------------------------------------------------------

antismash "$FASTA" \
    --genefinding-tool prodigal \
    --output-dir "$OUTDIR" \
    --output-basename "$ACC" \
    --cb-knownclusters \
    --clusterhmmer \
    --tigrfam \
    --asf \
    --rre \
    --pfam2go \
    --cpus 8 \
    > "$LOG" 2>&1


# ------------------------------------------------------------
# Completion check
# ------------------------------------------------------------

END=$(date +%s)

ELAPSED=$(( (END - START) / 60 ))

if [ -f "${OUTDIR}/${ACC}.json" ]; then

    echo "[$(date +'%Y-%m-%d %H:%M:%S')] DONE ${ACC} (${ELAPSED} min)" \
        | tee -a "$LOG_DIR/master.log"

else

    echo "[$(date +'%Y-%m-%d %H:%M:%S')] FAILED ${ACC} (${ELAPSED} min)" \
        | tee -a "$LOG_DIR/master.log"

    echo "Check:"
    echo "  $LOG"

    exit 1

fi
RUNNER_EOF

chmod +x "$RUNNER"


# ------------------------------------------------------------
# 7. Activate antiSMASH environment
# ------------------------------------------------------------

if ! command -v antismash >/dev/null 2>&1; then

    echo
    echo "ERROR: antiSMASH executable not found."
    echo
    echo "Activate the antiSMASH environment first, e.g.:"
    echo "  conda activate antismash"
    echo

    exit 1

fi


# ------------------------------------------------------------
# 8. Count genomes
# ------------------------------------------------------------

GENOME_COUNT=$(find "$GENOME_DIR" -maxdepth 1 -type f \
    \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) |
    wc -l)

echo
echo "Genome FASTA files detected: $GENOME_COUNT"

echo
echo "Entries in antiSMASH list:"
awk 'NF >= 2 && $1 !~ /^#/ {n++} END {print n+0}' "$LIST_FILE"


# ------------------------------------------------------------
# 9. Launch antiSMASH
# ------------------------------------------------------------

echo
echo "Launching antiSMASH with 2 concurrent genomes."
echo "Each antiSMASH process uses 8 CPUs."
echo

nohup bash -c '
    while IFS=$'\''\t'\'' read -r NAME ACC REST
    do
        [ -z "$NAME" ] && continue
        [[ "$NAME" =~ ^# ]] && continue

        "'"$RUNNER"'" "$NAME" "$ACC"

    done < "'"$LIST_FILE"'"
' > "$LOG_DIR/nohup_master.out" 2>&1 &

PID=$!

echo "$PID" > "$PROJECT_ROOT/run.pid"


# ------------------------------------------------------------
# 10. Final report
# ------------------------------------------------------------

echo
echo "============================================================"
echo "antiSMASH launched"
echo "============================================================"

echo
echo "PID:"
echo "  $PID"

echo
echo "PID file:"
echo "  $PROJECT_ROOT/run.pid"

echo
echo "Master log:"
echo "  $LOG_DIR/master.log"

echo
echo "nohup log:"
echo "  $LOG_DIR/nohup_master.out"

echo
echo "Output directory:"
echo "  $ANTISMASH_DIR"

echo
echo "Monitor with:"
echo "  tail -f $LOG_DIR/master.log"

echo
echo "Check process:"
echo "  ps -fp $PID"

echo
echo "============================================================"
