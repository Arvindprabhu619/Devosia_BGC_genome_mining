#!/usr/bin/env bash
# ============================================================
# scripts/04_run_checkm2.sh
# ============================================================
# Purpose: Estimate genome completeness and contamination.
#
# Inputs:   data/metadata/genome_manifest.tsv
# Outputs:
#   results/04_checkm2/checkm2_output/quality_report.tsv
#   results/04_checkm2/quality_report.tsv
#   results/04_checkm2/quality_passed.txt
#   results/04_checkm2/quality_excluded.txt
#
# Thresholds:
#   Standard: completeness >= 95%, contamination <= 5.0%
#   Borderline: completeness >= 95%, contamination 5.0-5.5%
#   Excluded: completeness < 95% OR contamination > 5.5%
#
# Runtime:  20-60 min
# ============================================================
set -euo pipefail
source "$(dirname "$0")/00_paths.sh"

LOG="${DIR_LOGS}/04_checkm2_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== [04] CheckM2: $(date) ==="
echo "CheckM2 DB: $CHECKM2_DB"

if [[ ! -f "$CHECKM2_DB" ]]; then
    echo "ERROR: CheckM2 DB not found at $CHECKM2_DB"
    exit 1
fi

# ---- Prepare input ----
INPUT_DIR="${RES_CHECKM2}/input_fastas"
mkdir -p "$INPUT_DIR"

echo "[04.1] Preparing input files..."
n=0
while IFS=$'\t' read -r acc fasta gff faa has_protein; do
    [[ "$acc" == "accession" ]] && continue
    [[ -z "$fasta" ]] && continue
    ln -sf "$fasta" "${INPUT_DIR}/${acc}.fna"
    n=$((n + 1))
done < "${DATA_META}/genome_manifest.tsv"
echo "  Linked $n FASTA files"

# ---- Run CheckM2 ----
if [[ ! -f "${RES_CHECKM2}/checkm2_output/quality_report.tsv" ]]; then
    echo "[04.2] Running CheckM2 (~20-60 min)..."
    $CHECKM2 predict \
        --input "$INPUT_DIR" \
        --output-directory "${RES_CHECKM2}/checkm2_output" \
        --database_path "$CHECKM2_DB" \
        --threads "$N_CPUS" \
        --force
else
    echo "[04.2] CheckM2 output already exists, skipping"
fi

# ---- Parse results ----
echo "[04.3] Parsing quality report..."

$PYTHON_BIN - <<PYEOF
import csv
from pathlib import Path

ROOT = Path("${RES_CHECKM2}")
report = ROOT / "checkm2_output" / "quality_report.tsv"

standard, borderline, excluded = [], [], []
rows_out = []

with open(report) as f:
    reader = csv.DictReader(f, delimiter="\\t")
    for row in reader:
        acc = row.get("Name", "").replace(".fna", "")
        try:
            comp = float(row.get("Completeness", 0))
            cont = float(row.get("Contamination", 0))
        except ValueError:
            continue

        if comp >= 95.0 and cont <= 5.0:
            status = "standard"
            standard.append(acc)
        elif comp >= 95.0 and cont <= 5.5:
            status = "borderline"
            borderline.append(acc)
        else:
            status = "excluded"
            excluded.append(acc)

        rows_out.append((acc, comp, cont, status))

with open(ROOT / "quality_report.tsv", "w") as f:
    f.write("accession\\tcompleteness\\tcontamination\\tstatus\\n")
    for r in rows_out:
        f.write(f"{r[0]}\\t{r[1]:.2f}\\t{r[2]:.2f}\\t{r[3]}\\n")

(ROOT / "quality_passed.txt").write_text("\\n".join(standard + borderline) + "\\n")
(ROOT / "quality_flags.tsv").write_text("\\n".join(borderline) + "\\n")
(ROOT / "quality_excluded.txt").write_text("\\n".join(excluded) + "\\n")

print(f"  Standard quality:  {len(standard)}")
print(f"  Borderline:        {len(borderline)}")
print(f"  Excluded:          {len(excluded)}")
print(f"  Total retained:    {len(standard) + len(borderline)}")
PYEOF

echo "=== [04] Done: $(date) ==="
