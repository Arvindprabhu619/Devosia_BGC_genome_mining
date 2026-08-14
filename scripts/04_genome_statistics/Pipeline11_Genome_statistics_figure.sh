#!/usr/bin/env bash

# ============================================================
# Pipeline 11 — Final genome statistics figures
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

FINAL="$PROJECT_ROOT/final_validation"
TABLE="$FINAL/genome_statistics/FINAL_110_Genome_Statistics.tsv"
FIG="$FINAL/figures"

mkdir -p "$FIG"

python - "$TABLE" "$FIG" <<'PY'
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

TABLE = Path(sys.argv[1])
FIG = Path(sys.argv[2])

df = pd.read_csv(TABLE, sep="\t")

print("=" * 70)
print("FINAL 110-GENOME ASSEMBLY FIGURES")
print("=" * 70)

print("Genomes:", len(df))

# ------------------------------------------------------------
# Figure 9 — Genome size distribution
# ------------------------------------------------------------

plt.figure(figsize=(8, 6))

plt.hist(
    df["Genome_Size_Mb"],
    bins=15
)

plt.xlabel("Genome size (Mb)")
plt.ylabel("Number of genomes")
plt.title("Genome size distribution of the 110 validated genomes")

plt.tight_layout()

plt.savefig(
    FIG / "Figure_9_Genome_size_distribution.png",
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# Figure 10 — GC content distribution
# ------------------------------------------------------------

plt.figure(figsize=(8, 6))

plt.hist(
    df["GC_Percent"],
    bins=15
)

plt.xlabel("GC content (%)")
plt.ylabel("Number of genomes")
plt.title("GC-content distribution of the 110 validated genomes")

plt.tight_layout()

plt.savefig(
    FIG / "Figure_10_GC_content_distribution.png",
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# Figure 11 — Genome size vs contig number
# ------------------------------------------------------------

plt.figure(figsize=(8, 6))

plt.scatter(
    df["Genome_Size_Mb"],
    df["Contigs"],
    alpha=0.8
)

plt.xlabel("Genome size (Mb)")
plt.ylabel("Number of contigs")
plt.title("Genome size versus assembly fragmentation")

plt.tight_layout()

plt.savefig(
    FIG / "Figure_11_Genome_size_vs_contigs.png",
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# Figure 12 — Genome size vs N50
# ------------------------------------------------------------

plt.figure(figsize=(8, 6))

plt.scatter(
    df["Genome_Size_Mb"],
    df["N50_bp"],
    alpha=0.8
)

plt.xlabel("Genome size (Mb)")
plt.ylabel("N50 (bp)")
plt.title("Genome size versus assembly N50")

plt.tight_layout()

plt.savefig(
    FIG / "Figure_12_Genome_size_vs_N50.png",
    dpi=300
)

plt.close()


print()
print("Figures generated:")

for f in sorted(FIG.glob("Figure_9*")):
    print(f.name)

for f in sorted(FIG.glob("Figure_10*")):
    print(f.name)

for f in sorted(FIG.glob("Figure_11*")):
    print(f.name)

for f in sorted(FIG.glob("Figure_12*")):
    print(f.name)

PY

echo
echo "============================================================"
echo "FIGURE FILES"
echo "============================================================"

ls -lh "$FIG"/Figure_{9,10,11,12}*.png
