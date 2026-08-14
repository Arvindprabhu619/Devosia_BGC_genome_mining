#!/usr/bin/env bash

# ============================================================
# Pipeline 10 — Final genome statistics
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

FINAL="$PROJECT_ROOT/final_validation"
GENOMES="$FINAL/final_genomes"
OUT="$FINAL/genome_statistics"

mkdir -p "$OUT"

echo "============================================================"
echo "FINAL 110 DEVOSIA GENOME STATISTICS"
echo "============================================================"

echo
echo "Genome directory:"
echo "$GENOMES"

echo
echo "FASTA files:"
find "$GENOMES" -maxdepth 1 -type f -name "*.fna" | wc -l

python - "$GENOMES" "$OUT" <<'PY'
import sys
from pathlib import Path
import csv
import statistics

GENOMES = Path(sys.argv[1])
OUT = Path(sys.argv[2])

OUT.mkdir(parents=True, exist_ok=True)

def fasta_stats(path):

    lengths = []
    gc = 0
    total = 0
    n_bases = 0
    current = 0

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current > 0:
                    lengths.append(current)
                current = 0
            else:
                seq = line.upper()
                current += len(seq)

                total += len(seq)
                gc += seq.count("G") + seq.count("C")
                n_bases += seq.count("N")

        if current > 0:
            lengths.append(current)

    lengths.sort(reverse=True)

    genome_size = sum(lengths)
    contigs = len(lengths)

    # N50
    half = genome_size / 2
    cumulative = 0
    n50 = 0
    l50 = 0

    for i, length in enumerate(lengths, 1):
        cumulative += length

        if cumulative >= half:
            n50 = length
            l50 = i
            break

    largest = lengths[0] if lengths else 0

    gc_percent = (
        (gc / (total - n_bases)) * 100
        if total > n_bases
        else 0
    )

    return {
        "Genome": path.name,
        "Genome_Size_bp": genome_size,
        "Genome_Size_Mb": genome_size / 1_000_000,
        "Contigs": contigs,
        "Largest_Contig_bp": largest,
        "N50_bp": n50,
        "L50": l50,
        "GC_Percent": gc_percent,
        "Ambiguous_N_bases": n_bases,
    }


files = sorted(GENOMES.glob("*.fna"))

rows = []

for i, path in enumerate(files, 1):

    stats = fasta_stats(path)
    rows.append(stats)

    print(
        f"[{i:3d}/{len(files)}] "
        f"{path.name} "
        f"{stats['Genome_Size_Mb']:.2f} Mb "
        f"{stats['Contigs']} contigs "
        f"N50={stats['N50_bp']:,}"
    )


outfile = OUT / "FINAL_110_Genome_Statistics.tsv"

fields = [
    "Genome",
    "Genome_Size_bp",
    "Genome_Size_Mb",
    "Contigs",
    "Largest_Contig_bp",
    "N50_bp",
    "L50",
    "GC_Percent",
    "Ambiguous_N_bases",
]

with open(outfile, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)


print()
print("============================================================")
print("SUMMARY")
print("============================================================")

sizes = [r["Genome_Size_Mb"] for r in rows]
contigs = [r["Contigs"] for r in rows]
n50s = [r["N50_bp"] for r in rows]
gc = [r["GC_Percent"] for r in rows]

print(f"Genomes analysed : {len(rows)}")

print(
    f"Genome size      : "
    f"{min(sizes):.2f}–{max(sizes):.2f} Mb"
)

print(
    f"Genome size mean : "
    f"{statistics.mean(sizes):.2f} Mb"
)

print(
    f"Contigs          : "
    f"{min(contigs)}–{max(contigs)}"
)

print(
    f"N50              : "
    f"{min(n50s):,}–{max(n50s):,} bp"
)

print(
    f"GC content       : "
    f"{min(gc):.2f}–{max(gc):.2f}%"
)

print()
print("Output:")
print(outfile)

PY

echo
echo "============================================================"
echo "OUTPUT"
echo "============================================================"

ls -lh "$OUT"

echo
echo "First 5 records:"
head -6 "$OUT/FINAL_110_Genome_Statistics.tsv"
