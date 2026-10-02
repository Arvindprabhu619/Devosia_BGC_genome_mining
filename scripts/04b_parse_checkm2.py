#!/usr/bin/env python3
"""Parse CheckM2 output into classification tables."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "03_checkm2"
REPORT = RES / "checkm2_output" / "quality_report.tsv"

standard, borderline, excluded = [], [], []
rows_out = []

with open(REPORT) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        name = row.get("Name", "")
        acc = name.replace(".fna", "")
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

with open(RES / "quality_report.tsv", "w") as f:
    f.write("accession\tcompleteness\tcontamination\tstatus\n")
    for r in rows_out:
        f.write(f"{r[0]}\t{r[1]:.2f}\t{r[2]:.2f}\t{r[3]}\n")

(RES / "quality_passed.txt").write_text("\n".join(standard + borderline) + "\n")
(RES / "quality_flags.tsv").write_text("\n".join(borderline) + "\n")
(RES / "quality_excluded.txt").write_text("\n".join(excluded) + "\n")

print(f"Standard quality:  {len(standard)}")
print(f"Borderline:        {len(borderline)}")
print(f"Excluded:          {len(excluded)}")
print(f"Total retained:    {len(standard) + len(borderline)}")
