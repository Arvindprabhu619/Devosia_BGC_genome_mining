#!/usr/bin/env bash
# ============================================================
# scripts/05_run_gtdbtk.sh
# ============================================================
# Purpose: Assign each quality-passing genome to GTDB R232 taxonomy.
#
# Inputs:   results/04_checkm2/quality_passed.txt
#           data/metadata/genome_manifest.tsv
# Outputs:
#   results/05_gtdbtk/gtdbtk_output/
#   results/05_gtdbtk/taxonomy_assignments.tsv
#   results/05_gtdbtk/final_genome_list.txt
#   results/05_gtdbtk/excluded_taxa.txt
#
# Tool:     GTDB-Tk 2.7.2
# Runtime:  20-60 min (with ANI pre-screening)
#
# Key fix:  Symlinks use 'query_' prefix to avoid ID collision with
#           GTDB reference genomes.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/05_gtdbtk_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [05] GTDB-Tk: $(date) ==="
echo "GTDB data: $GTDBTK_DATA"

if [[ ! -d "$GTDBTK_DATA" ]]; then
    echo "ERROR: GTDB data not found at $GTDBTK_DATA"
    exit 1
fi

# Ensure GTDBTK_DATA_PATH is exported (used internally by gtdbtk)
export GTDBTK_DATA_PATH="$GTDBTK_DATA"

# ---- Prepare input with 'query_' prefix ----
INPUT_DIR="${RES_GTDBTK}/input_genomes"
mkdir -p "$INPUT_DIR"

echo "[05.1] Linking quality-passed genomes with query_ prefix..."
n=0
while read -r acc; do
    [[ -z "$acc" ]] && continue
    fasta=$(awk -F'\t' -v a="$acc" '$1==a {print $2}' "${DATA_META}/genome_manifest.tsv")
    if [[ -n "$fasta" && -f "$fasta" ]]; then
        ln -sf "$fasta" "$INPUT_DIR/query_${acc}.fna"
        n=$((n + 1))
    fi
done < "${RES_CHECKM2}/quality_passed.txt"
echo "  Linked $n genomes"

# ---- Run classify_wf ----
if [[ ! -f "${RES_GTDBTK}/gtdbtk_output/classify/gtdbtk.bac120.summary.tsv" ]]; then
    echo "[05.2] Running GTDB-Tk classify_wf (20-60 min)..."
    $GTDBTK classify_wf \
        --genome_dir "$INPUT_DIR" \
        --out_dir "${RES_GTDBTK}/gtdbtk_output" \
        --extension fna \
        --cpus "$N_CPUS" \
        --pplacer_cpus "$N_CPUS"
else
    echo "[05.2] GTDB-Tk output already exists, skipping"
fi

# ---- Parse taxonomy ----
echo "[05.3] Parsing taxonomy..."

$PYTHON_BIN - <<PYEOF
import csv
from pathlib import Path
from collections import Counter

RES = Path("${RES_GTDBTK}")
summary = RES / "gtdbtk_output" / "classify" / "gtdbtk.bac120.summary.tsv"

rows = []
with open(summary) as f:
    reader = csv.DictReader(f, delimiter="\\t")
    for row in reader:
        rows.append(row)

print(f"Total taxonomy rows: {len(rows)}")

def get_genus(tax):
    for part in tax.split(";"):
        if part.startswith("g__"):
            return part.replace("g__", "")
    return "Unknown"

devosia, devosia_a, other = [], [], []
genus_counts = Counter()

for row in rows:
    acc = row["user_genome"].replace(".fna", "").replace("query_", "")
    tax = row["classification"]
    genus = get_genus(tax)
    genus_counts[genus] += 1

    if genus == "Devosia_A":
        devosia_a.append(acc)
    elif genus == "Devosia":
        devosia.append(acc)
    else:
        other.append((acc, genus))

print("\\n=== Genus distribution ===")
for g, c in genus_counts.most_common():
    print(f"  {g:30s} {c}")

with open(RES / "taxonomy_assignments.tsv", "w") as f:
    f.write("accession\\ttaxonomy\\tgenus\\tstatus\\n")
    for row in rows:
        acc = row["user_genome"].replace(".fna", "").replace("query_", "")
        tax = row["classification"]
        genus = get_genus(tax)
        status = "Devosia" if genus == "Devosia" else ("Devosia_A" if genus == "Devosia_A" else "Excluded")
        f.write(f"{acc}\\t{tax}\\t{genus}\\t{status}\\n")

with open(RES / "final_genome_list.txt", "w") as f:
    for acc in devosia + devosia_a:
        f.write(acc + "\\n")

with open(RES / "excluded_taxa.txt", "w") as f:
    for acc, genus in other:
        f.write(f"{acc}\\t{genus}\\n")

print(f"\\n=== Final ===")
print(f"  Devosia:      {len(devosia)}")
print(f"  Devosia_A:    {len(devosia_a)}")
print(f"  Excluded:     {len(other)}")
print(f"  FINAL:        {len(devosia) + len(devosia_a)}")
PYEOF

echo "=== [05] Done: $(date) ==="
