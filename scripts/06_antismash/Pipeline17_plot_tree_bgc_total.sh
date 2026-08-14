#!/usr/bin/env bash

# ============================================================
# Pipeline 17 — Phylogeny + total BGC abundance
#
# Purpose:
#   Map total predicted BGC abundance onto the final rooted
#   110-genome Devosia phylogeny.
#
# Inputs:
#   phylogeny_110/Devosia_110_final_rooted_pruned.tree
#   bgc_counts_per_genome_fixed.tsv
#
# Output:
#   manuscript_outputs/Figure5_tree_BGC_total.png
#
# Portability:
#   Repository root is inferred from this script's location.
#   No hard-coded HOME or project-specific paths.
# ============================================================

set -euo pipefail

## -----------------------------------------------------------
## 1. Determine repository root
## -----------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

## -----------------------------------------------------------
## 2. Project paths
## -----------------------------------------------------------

PHYLO_DIR="$PROJECT_ROOT/phylogeny_110"
TREE_FILE="$PHYLO_DIR/Devosia_110_final_rooted_pruned.tree"
BGC_FILE="$PROJECT_ROOT/bgc_counts_per_genome_fixed.tsv"
OUTDIR="$PROJECT_ROOT/manuscript_outputs"
PYTHON_SCRIPT="$PHYLO_DIR/plot_tree_bgc_total.py"
OUTPUT="$OUTDIR/Figure5_tree_BGC_total.png"

## -----------------------------------------------------------
## 3. Input validation
## -----------------------------------------------------------

if [ ! -f "$TREE_FILE" ]; then
    echo "ERROR: Final rooted phylogeny not found:"
    echo "  $TREE_FILE"
    exit 1
fi

if [ ! -f "$BGC_FILE" ]; then
    echo "ERROR: Per-genome BGC count file not found:"
    echo "  $BGC_FILE"
    exit 1
fi

mkdir -p "$OUTDIR"

## -----------------------------------------------------------
## 4. Report configuration
## -----------------------------------------------------------

echo "============================================================"
echo "Pipeline 17 — Phylogeny + total BGC abundance"
echo "============================================================"

echo
echo "Project root:"
echo "  $PROJECT_ROOT"

echo
echo "Phylogeny:"
echo "  $TREE_FILE"

echo
echo "BGC counts:"
echo "  $BGC_FILE"

echo
echo "Output directory:"
echo "  $OUTDIR"

echo

TREE_FILE="$TREE_FILE" \
BGC_FILE="$BGC_FILE" \
OUTPUT="$OUTPUT" \
python3 - <<'PY'

from pathlib import Path
import os
import sys

import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from Bio import Phylo


# ============================================================
# Paths supplied by shell wrapper
# ============================================================

tree_file = Path(os.environ["TREE_FILE"])
bgc_file = Path(os.environ["BGC_FILE"])
output = Path(os.environ["OUTPUT"])

output.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. Read final rooted phylogeny
# ============================================================

print("Reading phylogeny...")

tree = Phylo.read(str(tree_file), "newick")

terminals = tree.get_terminals()

print(f"Tree terminals: {len(terminals)}")

EXPECTED_TIPS = 110

if len(terminals) != EXPECTED_TIPS:
    raise SystemExit(
        f"ERROR: Expected {EXPECTED_TIPS} tree tips, "
        f"found {len(terminals)}."
    )


# ============================================================
# 2. Read per-genome BGC abundance
# ============================================================

print("Reading per-genome BGC counts...")

bgc = pd.read_csv(
    bgc_file,
    sep="\t",
    names=["accession", "n_bgcs"],
    header=None
)

if bgc.empty:
    raise SystemExit("ERROR: BGC count table is empty.")

bgc["accession"] = bgc["accession"].astype(str)
bgc["n_bgcs"] = pd.to_numeric(
    bgc["n_bgcs"],
    errors="coerce"
)

if bgc["n_bgcs"].isna().any():
    raise SystemExit(
        "ERROR: Non-numeric BGC counts detected."
    )

print(f"BGC table rows: {len(bgc)}")

if bgc["accession"].duplicated().any():
    duplicated = bgc.loc[
        bgc["accession"].duplicated(),
        "accession"
    ].tolist()

    raise SystemExit(
        "ERROR: Duplicate genome accessions detected:\n"
        + "\n".join(duplicated)
    )


# ============================================================
# 3. Match tree tips to BGC counts
# ============================================================

import re

def normalize_accession(x):
    """
    Normalize genome identifiers to the stable GCF accession.

    Tree example:
        GCF_000970435.1_ASM97043v1_genomic

    BGC table example:
        GCF_000970435.1

    Both normalize to:
        GCF_000970435.1
    """
    x = str(x)
    m = re.match(r"(GCF_\d+\.\d+)", x)
    return m.group(1) if m else x


tree_accessions = [str(t.name) for t in terminals]

tree_norm = {
    normalize_accession(acc): acc
    for acc in tree_accessions
}

bgc_norm = bgc.copy()
bgc_norm["normalized_accession"] = bgc_norm["accession"].map(
    normalize_accession
)

if bgc_norm["normalized_accession"].duplicated().any():
    duplicated = bgc_norm.loc[
        bgc_norm["normalized_accession"].duplicated(),
        "normalized_accession"
    ].tolist()

    raise SystemExit(
        "ERROR: Duplicate normalized genome accessions detected:\\n"
        + "\\n".join(duplicated)
    )

acc_to_n = dict(
    zip(
        bgc_norm["normalized_accession"],
        bgc_norm["n_bgcs"]
    )
)

missing = [
    acc for acc in tree_accessions
    if normalize_accession(acc) not in acc_to_n
]

extra = [
    acc for acc in bgc["accession"]
    if normalize_accession(acc) not in tree_norm
]

print(
    f"Normalized tree accessions: {len(tree_norm)}"
)

print(
    f"Normalized BGC accessions:  "
    f"{len(bgc_norm['normalized_accession'].unique())}"
)

print(
    f"Successfully matched:      "
    f"{len(tree_accessions) - len(missing)}/{len(tree_accessions)}"
)

print(f"Tree tips missing BGC counts: {len(missing)}")
print(f"BGC accessions absent from tree: {len(extra)}")

if missing:
    print("\nMissing BGC counts:")
    for acc in missing:
        print(f"  {acc}")

    raise SystemExit(
        "ERROR: Not all tree genomes have BGC counts."
    )

if len(bgc) != EXPECTED_TIPS:
    raise SystemExit(
        f"ERROR: Expected {EXPECTED_TIPS} BGC-count rows, "
        f"found {len(bgc)}."
    )


# ============================================================
# 4. Assign tree coordinates
# ============================================================

print("Calculating tree coordinates...")

y_pos = {
    terminal: i + 1
    for i, terminal in enumerate(terminals)
}


def assign_positions(clade, parent_x=0.0):

    x = parent_x + (clade.branch_length or 0.0)

    clade._x = x

    if clade.is_terminal():

        clade._y = y_pos[clade]

    else:

        for child in clade.clades:
            assign_positions(child, x)

        ys = [
            child._y
            for child in clade.clades
        ]

        clade._y = sum(ys) / len(ys)

    return x


assign_positions(tree.root)


# ============================================================
# 5. Draw phylogeny
# ============================================================

def draw_clade(clade, ax):

    for child in clade.clades:

        ax.plot(
            [clade._x, child._x],
            [child._y, child._y],
            linewidth=0.7
        )

        draw_clade(child, ax)

    if not clade.is_terminal():

        ys = [
            child._y
            for child in clade.clades
        ]

        ax.plot(
            [clade._x, clade._x],
            [min(ys), max(ys)],
            linewidth=0.7
        )


# ============================================================
# 6. Create figure
# ============================================================

fig, (ax_tree, ax_bar) = plt.subplots(
    1,
    2,
    figsize=(16, 18),
    gridspec_kw={
        "width_ratios": [3.2, 1.0],
        "wspace": 0.05
    }
)


# ============================================================
# 7. Tree panel
# ============================================================

draw_clade(tree.root, ax_tree)

ax_tree.set_title(
    "Devosia phylogeny (n = 110)",
    fontsize=14,
    fontweight="bold",
    pad=15
)

ax_tree.set_xlabel(
    "Evolutionary distance (substitutions/site)",
    fontsize=11
)

ax_tree.set_ylabel("Genomes", fontsize=11)

ax_tree.set_yticks([])

ax_tree.spines["top"].set_visible(False)
ax_tree.spines["right"].set_visible(False)


# ============================================================
# 8. BGC abundance panel
# ============================================================

bar_values = [
    acc_to_n[normalize_accession(terminal.name)]
    for terminal in terminals
]

y_values = [
    terminal._y
    for terminal in terminals
]

ax_bar.barh(
    y_values,
    bar_values,
    height=0.72
)

ax_bar.set_xlabel(
    "Total BGCs",
    fontsize=11
)

ax_bar.set_title(
    "BGC abundance",
    fontsize=12,
    fontweight="bold"
)

ax_bar.set_yticks([])

ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)


# ============================================================
# 9. Axis alignment
# ============================================================

ax_tree.set_ylim(
    ax_bar.get_ylim()
)

ax_tree.set_xlim(
    left=0
)

ax_bar.set_ylim(
    ax_tree.get_ylim()
)


# ============================================================
# 10. Add summary statistics
# ============================================================

mean_bgc = bgc["n_bgcs"].mean()
sd_bgc = bgc["n_bgcs"].std()
min_bgc = bgc["n_bgcs"].min()
max_bgc = bgc["n_bgcs"].max()

fig.suptitle(
    (
        "Phylogenetic distribution of total biosynthetic gene cluster abundance\n"
        f"Mean = {mean_bgc:.2f} ± {sd_bgc:.2f}; "
        f"range = {int(min_bgc)}–{int(max_bgc)} BGCs/genome"
    ),
    fontsize=15,
    fontweight="bold",
    y=0.995
)


# ============================================================
# 11. Save publication figure
# ============================================================

plt.savefig(
    output,
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print()
print("============================================================")
print("Pipeline 17 Python analysis completed")
print("============================================================")
print()
print(f"Tree genomes:       {len(terminals)}")
print(f"BGC-count genomes:  {len(bgc)}")
print(f"Mean BGC count:     {mean_bgc:.3f}")
print(f"SD:                 {sd_bgc:.3f}")
print(f"Minimum:            {int(min_bgc)}")
print(f"Maximum:            {int(max_bgc)}")
print()
print(f"Output:")
print(f"  {output}")
print()


# ============================================================
# 12. Validate output
# ============================================================

if not output.exists():
    raise SystemExit(
        "ERROR: Expected figure was not created."
    )

if output.stat().st_size == 0:
    raise SystemExit(
        "ERROR: Figure file is empty."
    )

print(
    f"Output size: {output.stat().st_size / 1024:.1f} KB"
)

PY

echo
echo "============================================================"
echo "Pipeline 17 completed successfully."
echo "============================================================"

echo
echo "Figure:"
echo "  $OUTPUT"

echo
echo "Output directory:"
echo "  $OUTDIR"

echo
ls -lh "$OUTPUT"

