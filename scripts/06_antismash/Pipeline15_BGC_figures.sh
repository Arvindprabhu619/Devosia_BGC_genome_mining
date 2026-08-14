#!/usr/bin/env bash

# ============================================================
# Pipeline 15 — BGC summary tables and figures
#
# Purpose:
#   Generate manuscript-level BGC class summaries, genome-by-
#   class matrix, and BGC distribution figures from the
#   antiSMASH-derived region table.
#
# Inputs:
#   bgc_regions_detailed.tsv
#   bgc_counts_per_genome_fixed.tsv
#
# Outputs:
#   manuscript_outputs/
#
# Portability:
#   Repository root is inferred from this script's location.
#   No hard-coded HOME or project-specific paths.
#
# Dependencies:
#   Python 3
#   pandas
#   matplotlib
# ============================================================

set -euo pipefail

## ---- Determine repository root ----

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

## ---- Input/output paths ----

REGIONS="$PROJECT_ROOT/bgc_regions_detailed.tsv"
PER_GENOME="$PROJECT_ROOT/bgc_counts_per_genome_fixed.tsv"
OUTDIR="$PROJECT_ROOT/manuscript_outputs"

## ---- Input validation ----

if [ ! -f "$REGIONS" ]; then
    echo "ERROR: BGC region table not found:"
    echo "  $REGIONS"
    echo
    echo "Run Pipeline 14 first."
    exit 1
fi

if [ ! -f "$PER_GENOME" ]; then
    echo "ERROR: Per-genome BGC count table not found:"
    echo "  $PER_GENOME"
    exit 1
fi

## ---- Python dependency check ----

python3 - <<'PY'
try:
    import pandas
    import matplotlib
except ImportError as e:
    raise SystemExit(
        "ERROR: Required Python package is missing: "
        + str(e)
        + "\nInstall pandas and matplotlib in the analysis environment."
    )
PY

mkdir -p "$OUTDIR"

## ---- Report configuration ----

echo "============================================================"
echo "Pipeline 15 — BGC manuscript outputs"
echo "============================================================"

echo
echo "Project root:"
echo "  $PROJECT_ROOT"

echo
echo "BGC region table:"
echo "  $REGIONS"

echo
echo "Per-genome BGC table:"
echo "  $PER_GENOME"

echo
echo "Output directory:"
echo "  $OUTDIR"

## ---- Generate manuscript outputs ----

python3 - "$REGIONS" "$PER_GENOME" "$OUTDIR" <<'PY'
#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


REGIONS_FILE = Path(sys.argv[1])
PER_GENOME_FILE = Path(sys.argv[2])
OUTDIR = Path(sys.argv[3])

OUTDIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. READ INPUT TABLES
# ============================================================

regions = pd.read_csv(
    REGIONS_FILE,
    sep="\t"
)

per_genome = pd.read_csv(
    PER_GENOME_FILE,
    sep="\t",
    names=["accession", "n_bgcs"],
    header=None
)


required_region_columns = {
    "accession",
    "contig",
    "products"
}

missing_region_columns = (
    required_region_columns
    - set(regions.columns)
)

if missing_region_columns:
    raise SystemExit(
        "ERROR: Missing columns from BGC region table: "
        + ", ".join(sorted(missing_region_columns))
    )


required_per_genome_columns = {
    "accession",
    "n_bgcs"
}

if set(per_genome.columns) != required_per_genome_columns:
    raise SystemExit(
        "ERROR: Unexpected per-genome BGC table structure."
    )


# ============================================================
# 2. EXPAND BGC PRODUCT CLASSES
# ============================================================

regions["products"] = (
    regions["products"]
    .fillna("unknown")
)

regions_exploded = (
    regions.assign(
        product=regions["products"].str.split(";")
    )
    .explode("product")
)

n_regions_total = len(regions)
n_tag_instances = len(regions_exploded)

print(
    f"Total distinct BGC regions: {n_regions_total}"
)

print(
    f"Total class-tag instances: {n_tag_instances} "
    f"({n_tag_instances - n_regions_total} "
    f"from hybrid/multi-class regions)"
)


# ============================================================
# 3. BGC CLASS SUMMARY
# ============================================================

class_counts = (
    regions_exploded["product"]
    .value_counts()
    .rename_axis("BGC class")
    .reset_index(name="Count")
)

class_counts["% of total regions"] = (
    100
    * class_counts["Count"]
    / n_regions_total
).round(1)

class_counts.to_csv(
    OUTDIR / "Table1_BGC_class_summary.tsv",
    sep="\t",
    index=False
)

print()
print(
    f"Wrote {OUTDIR / 'Table1_BGC_class_summary.tsv'}"
)

print(
    class_counts.to_string(index=False)
)


# ============================================================
# 4. GENOME × BGC CLASS MATRIX
# ============================================================

matrix = (
    regions_exploded
    .groupby(["accession", "product"])
    .size()
    .unstack(fill_value=0)
)

matrix = matrix.reindex(
    per_genome["accession"],
    fill_value=0
)

matrix.to_csv(
    OUTDIR / "TableS1_genome_by_class_matrix.tsv",
    sep="\t"
)

print()
print(
    f"Wrote {OUTDIR / 'TableS1_genome_by_class_matrix.tsv'} "
    f"({matrix.shape[0]} genomes x "
    f"{matrix.shape[1]} classes)"
)


# ============================================================
# 5. FIGURE 1 — BGC CLASS DISTRIBUTION
# ============================================================

cc_sorted = class_counts.sort_values(
    "Count",
    ascending=True
)

fig, ax = plt.subplots(
    figsize=(9, 6)
)

colors = plt.cm.viridis_r(
    np.linspace(
        0,
        1,
        len(cc_sorted)
    )
)

ax.barh(
    cc_sorted["BGC class"],
    cc_sorted["Count"],
    color=colors
)

ax.set_xlabel(
    "Number of BGC regions",
    fontsize=11
)

ax.set_title(
    f"Distribution of BGC classes across "
    f"{matrix.shape[0]} Devosia genomes "
    f"(n = {n_regions_total} regions)",
    fontsize=12
)

ax.spines[
    ["top", "right"]
].set_visible(False)

for i, v in enumerate(
    cc_sorted["Count"]
):
    ax.text(
        v + 1.5,
        i,
        str(v),
        va="center",
        fontsize=8
    )

plt.tight_layout()

figure1 = (
    OUTDIR /
    "Figure1_BGC_class_distribution.png"
)

plt.savefig(
    figure1,
    dpi=300
)

plt.close()

print(
    f"Wrote {figure1}"
)


# ============================================================
# 6. FIGURE 2 — PER-GENOME BGC HISTOGRAM
# ============================================================

counts = per_genome["n_bgcs"]

mean = counts.mean()
sd = counts.std()

bins = range(
    int(counts.min()),
    int(counts.max()) + 2
)

fig, ax = plt.subplots(
    figsize=(7, 5)
)

ax.hist(
    counts,
    bins=bins,
    align="left",
    color="#4C72B0",
    edgecolor="black",
    rwidth=0.7
)

ax.set_xlabel(
    "BGCs per genome",
    fontsize=11
)

ax.set_ylabel(
    "Number of genomes",
    fontsize=11
)

ax.set_title(
    f"Per-genome BGC count distribution "
    f"(n = {len(counts)} genomes)\n"
    f"Mean = {mean:.2f} ± {sd:.2f} SD",
    fontsize=12
)

ax.set_xticks(
    list(bins)[:-1]
)

ax.spines[
    ["top", "right"]
].set_visible(False)

plt.tight_layout()

figure2 = (
    OUTDIR /
    "Figure2_per_genome_BGC_histogram.png"
)

plt.savefig(
    figure2,
    dpi=300
)

plt.close()

print(
    f"Wrote {figure2}"
)


# ============================================================
# 7. FIGURE 3 — GENOME × BGC CLASS HEATMAP
# ============================================================

fig, ax = plt.subplots(
    figsize=(
        12,
        max(
            6,
            matrix.shape[0] * 0.08
        )
    )
)

im = ax.imshow(
    matrix.values,
    aspect="auto",
    cmap="viridis"
)

ax.set_xticks(
    range(matrix.shape[1])
)

ax.set_xticklabels(
    matrix.columns,
    rotation=90,
    fontsize=7
)

ax.set_yticks([])

ax.set_ylabel(
    f"Genomes (n = {matrix.shape[0]})",
    fontsize=11
)

ax.set_title(
    "BGC class presence/count per genome",
    fontsize=12
)

cbar = plt.colorbar(
    im,
    ax=ax,
    shrink=0.6
)

cbar.set_label(
    "Count per genome",
    fontsize=9
)

plt.tight_layout()

figure3 = (
    OUTDIR /
    "Figure3_genome_by_class_heatmap.png"
)

plt.savefig(
    figure3,
    dpi=300
)

plt.close()

print(
    f"Wrote {figure3}"
)


# ============================================================
# 8. FINAL SUMMARY
# ============================================================

print()
print(
    "============================================================"
)

print(
    "Pipeline 15 Python analysis completed."
)

print(
    "============================================================"
)

print(
    f"Distinct BGC regions : {n_regions_total}"
)

print(
    f"Class-tag instances  : {n_tag_instances}"
)

print(
    f"Genomes              : {matrix.shape[0]}"
)

print(
    f"BGC classes          : {matrix.shape[1]}"
)

print()
print(
    f"Output directory: {OUTDIR}"
)

PY

echo
echo "============================================================"
echo "Pipeline 15 completed."
echo "============================================================"

echo
echo "Generated outputs:"
find "$OUTDIR" -maxdepth 1 -type f -printf '  %f\n' | sort
