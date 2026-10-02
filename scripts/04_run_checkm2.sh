#!/usr/bin/env bash
# ============================================================
# scripts/04_run_checkm2.sh
#
# Purpose: Run CheckM2 on all curated genomes to estimate
#          completeness and contamination.
#
# Inputs:  data/metadata/genome_manifest.tsv
# Outputs: results/03_checkm2/quality_report.tsv
#          results/03_checkm2/quality_passed.txt
#          results/03_checkm2/quality_flags.tsv (borderline)
#          results/03_checkm2/quality_excluded.txt
#
# Thresholds:
#   Standard quality:  completeness >= 95%, contamination <= 5.0%
#   Borderline:        completeness >= 95%, contamination 5.0-5.5%
#   Excluded:          completeness < 95% OR contamination > 5.5%
#
# Runtime: 1-3 hours
# ============================================================

set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/04_checkm2_$(date +%Y%m%d_%H%M).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [04] CheckM2: $(date) ==="
echo "CheckM2 DB: $CHECKM2_DB"

# Verify CheckM2 DB exists
if [[ ! -f "$CHECKM2_DB" ]]; then
    echo "ERROR: CheckM2 DB not found at $CHECKM2_DB"
    echo "Wait for the database download to finish."
    exit 1
fi

mkdir -p "$RES_CHECKM2"

# Build input directory: symlink all FASTA files
INPUT_DIR="${RES_CHECKM2}/input_fastas"
mkdir -p "$INPUT_DIR"

echo "[04.1] Preparing input directory..."
n=0
while IFS=$'\t' read -r acc fasta gff faa has_protein; do
    [[ "$acc" == "accession" ]] && continue
    [[ -z "$fasta" ]] && continue
    ln -sf "$fasta" "$INPUT_DIR/${acc}.fna"
    n=$((n+1))
done < "${DATA_META}/genome_manifest.tsv"
echo "  Linked $n FASTA files"

# Run CheckM2
echo "[04.2] Running CheckM2 (this takes 1-3 hours)..."
$CHECKM2 predict \
    --input "$INPUT_DIR" \
    --output-directory "${RES_CHECKM2}/checkm2_output" \
    --database_path "$CHECKM2_DB" \
    --threads "$N_CPUS" \
    --force

echo "[04.3] Parsing quality report..."

# Parse quality report
python3 - <<'PYEOF' || ~/miniconda3/envs/checkm2/bin/python - <<'PYEOF2'
import csv
from pathlib import Path

root = Path("${RES_CHECKM2}")
report = root / "checkm2_output" / "quality_report.tsv"

# Read CheckM2 output
standard, borderline, excluded = [], [], []
rows_out = []

with open(report) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        name = row.get("Name") or row.get("name") or ""
        # Strip .fna extension
        acc = name.replace(".fna", "")
        try:
            completeness = float(row.get("Completeness") or row.get("completeness") or 0)
            contamination = float(row.get("Contamination") or row.get("contamination") or 0)
        except ValueError:
            continue

        if completeness >= 95.0 and contamination <= 5.0:
            status = "standard"
            standard.append(acc)
        elif completeness >= 95.0 and contamination <= 5.5:
            status = "borderline"
            borderline.append(acc)
        else:
            status = "excluded"
            excluded.append(acc)

        rows_out.append((acc, completeness, contamination, status))

# Write outputs
with open(root / "quality_report.tsv", "w") as f:
    f.write("accession\tcompleteness\tcontamination\tstatus\n")
    for r in rows_out:
        f.write(f"{r[0]}\t{r[1]:.2f}\t{r[2]:.2f}\t{r[3]}\n")

(root / "quality_passed.txt").write_text("\n".join(standard + borderline) + "\n")
(root / "quality_flags.tsv").write_text("\n".join(borderline) + "\n")
(root / "quality_excluded.txt").write_text("\n".join(excluded) + "\n")

print(f"Standard quality:  {len(standard)}")
print(f"Borderline:        {len(borderline)}")
print(f"Excluded:          {len(excluded)}")
print(f"Total retained:    {len(standard) + len(borderline)}")
PYEOF

echo "=== [04] Done: $(date) ==="
