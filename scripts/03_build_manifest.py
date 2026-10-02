#!/usr/bin/env python3
"""
03_build_manifest.py

Builds a manifest of all downloaded genomes:
  - Accession
  - FASTA path
  - GFF3 path (may be missing)
  - FAA path (may be missing)
  - Has_protein flag
Outputs:
  - data/metadata/genome_manifest.tsv
  - data/metadata/genomes_with_protein.txt
  - data/metadata/genomes_without_protein.txt
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw_genomes" / "ncbi_dataset" / "data"
META = ROOT / "data" / "metadata"

rows = []
missing_protein = []
with_protein = []

for acc_dir in sorted(RAW.iterdir()):
    if not acc_dir.is_dir():
        continue
    acc = acc_dir.name

    # Find genomic FASTA (the assembly, not cds_from_genomic)
    fasta_files = [f for f in acc_dir.glob("*.fna") if "cds_from_genomic" not in f.name]
    gff_files = list(acc_dir.glob("*.gff"))
    faa_files = list(acc_dir.glob("*.faa"))

    fasta = fasta_files[0] if fasta_files else None
    gff = gff_files[0] if gff_files else None
    faa = faa_files[0] if faa_files else None

    has_protein = gff is not None and faa is not None

    rows.append({
        "accession": acc,
        "fasta": str(fasta) if fasta else "",
        "gff": str(gff) if gff else "",
        "faa": str(faa) if faa else "",
        "has_protein": "yes" if has_protein else "no",
    })

    if has_protein:
        with_protein.append(acc)
    else:
        missing_protein.append(acc)

# Write manifest
header = "accession\tfasta\tgff\tfaa\thas_protein"
lines = [header]
for r in rows:
    lines.append(f"{r['accession']}\t{r['fasta']}\t{r['gff']}\t{r['faa']}\t{r['has_protein']}")
(META / "genome_manifest.tsv").write_text("\n".join(lines))

# Write per-category accession lists
(META / "genomes_with_protein.txt").write_text("\n".join(with_protein) + "\n")
(META / "genomes_without_protein.txt").write_text("\n".join(missing_protein) + "\n")

print(f"Total genomes in manifest:  {len(rows)}")
print(f"With protein/GFF3:          {len(with_protein)}")
print(f"Without protein/GFF3:       {len(missing_protein)}")
print()
print(f"Manifest:                   {META / 'genome_manifest.tsv'}")
print(f"With protein list:          {META / 'genomes_with_protein.txt'}")
print(f"Without protein list:       {META / 'genomes_without_protein.txt'}")
