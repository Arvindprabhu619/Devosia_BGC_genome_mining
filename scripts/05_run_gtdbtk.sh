#!/usr/bin/env bash
# ============================================================
# scripts/05_run_gtdbtk.sh
# Purpose: Assign each quality-passing genome to GTDB R232 taxonomy.
# Inputs:  results/03_checkm2/quality_passed.txt
# Outputs: results/04_gtdbtk/gtdbtk_output/
#          results/04_gtdbtk/taxonomy_assignments.tsv
#          results/04_gtdbtk/final_genome_list.txt
#          results/04_gtdbtk/excluded_taxa.txt
# Tool: GTDB-Tk 2.7.2
# Runtime: 4-6 hours
# ============================================================

set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/05_gtdbtk_$(date +%Y%m%d_%H%M).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [05] GTDB-Tk: $(date) ==="
echo "GTDB data: $GTDBTK_DATA"

if [[ ! -d "$GTDBTK_DATA" ]]; then
    echo "ERROR: GTDB data not found at $GTDBTK_DATA"
    exit 1
fi

# GTDB-Tk needs GTDBTK_DATA_PATH
export GTDBTK_DATA_PATH="$GTDBTK_DATA"

mkdir -p "$RES_GTDBTK"

# Prepare input directory
INPUT_DIR="${RES_GTDBTK}/input_genomes"
mkdir -p "$INPUT_DIR"

echo "[05.1] Linking quality-passed genomes..."
n=0
while read -r acc; do
    [[ -z "$acc" ]] && continue
    fasta=$(awk -F'\t' -v a="$acc" '$1==a {print $2}' "${DATA_META}/genome_manifest.tsv")
    if [[ -n "$fasta" && -f "$fasta" ]]; then
        ln -sf "$fasta" "$INPUT_DIR/${acc}.fna"
        n=$((n+1))
    fi
done < "${RES_CHECKM2}/quality_passed.txt"
echo "  Linked $n genomes"

# Run GTDB-Tk classify_wf
echo "[05.2] Running GTDB-Tk classify_wf (4-6 hours)..."

$GTDBTK classify_wf \
    --genome_dir "$INPUT_DIR" \
    --out_dir "${RES_GTDBTK}/gtdbtk_output" \
    --extension fna \
    --cpus "$N_CPUS" \
    --pplacer_cpus "$N_CPUS"

echo "[05.3] Parsing taxonomy assignments..."

# Parse GTDB-Tk summary
python3 - <<'PYEOF' || ~/miniconda3/envs/checkm2/bin/python - <<'PYEOF2'
from pathlib import Path

gtdbtk_dir = Path("${RES_GTDBTK}/gtdbtk_output")
summary_file = gtdbtk_dir / "classify" / "gtdbtk.bac120.summary.tsv"

devosia, devosia_a, excluded = [], [], []
rows = []

with open(summary_file) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        acc = parts[0].replace(".fna", "")
        taxonomy = parts[1]
        rows.append((acc, taxonomy))

        if "Devosia_A" in taxonomy:
            devosia_a.append(acc)
        elif "Devosia" in taxonomy:
            devosia.append(acc)
        else:
            excluded.append(acc)

with open("${RES_GTDBTK}/taxonomy_assignments.tsv", "w") as f:
    f.write("accession\ttaxonomy\tstatus\n")
    for acc, tax in rows:
        if acc in devosia:
            status = "Devosia"
        elif acc in devosia_a:
            status = "Devosia_A"
        else:
            status = "Excluded"
        f.write(f"{acc}\t{tax}\t{status}\n")

with open("${RES_GTDBTK}/final_genome_list.txt", "w") as f:
    for acc in devosia + devosia_a:
        f.write(acc + "\n")

with open("${RES_GTDBTK}/excluded_taxa.txt", "w") as f:
    for acc in excluded:
        f.write(acc + "\n")

print(f"Devosia:       {len(devosia)}")
print(f"Devosia_A:     {len(devosia_a)}")
print(f"Excluded:      {len(excluded)}")
print(f"Final list:    {len(devosia) + len(devosia_a)}")
PYEOF

echo "=== [05] Done: $(date) ==="
