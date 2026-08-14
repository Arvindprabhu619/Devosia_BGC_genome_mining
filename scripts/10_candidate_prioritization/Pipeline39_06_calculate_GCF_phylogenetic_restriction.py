#!/usr/bin/env python3

import re
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Phylo

TREE = Path("phylogeny_110/Devosia_110_final_rooted_pruned.tree")
EVIDENCE = "BGC_candidate_evidence_c0.7_corrected.tsv"

OUT_GCF = "GCF_phylogenetic_restriction_c0.7.tsv"
OUT_BGC = "BGC_candidate_evidence_phylogeny_c0.7.tsv"


def extract_gcf(x):
    m = re.search(r"(GCF_\d+\.\d+)", str(x))
    return m.group(1) if m else None


def mean_pairwise_distance(tree, names):
    names = list(names)

    if len(names) < 2:
        return 0.0

    vals = []

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            try:
                d = tree.distance(names[i], names[j])
                vals.append(float(d))
            except Exception:
                pass

    return float(np.mean(vals)) if vals else np.nan


def median_pairwise_distance(tree, names):
    names = list(names)

    if len(names) < 2:
        return 0.0

    vals = []

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            try:
                d = tree.distance(names[i], names[j])
                vals.append(float(d))
            except Exception:
                pass

    return float(np.median(vals)) if vals else np.nan


print("=" * 100)
print("GCF PHYLOGENETIC RESTRICTION ANALYSIS")
print("=" * 100)

# ============================================================
# TREE
# ============================================================

if not TREE.exists():
    raise SystemExit(f"ERROR: missing tree: {TREE}")

tree = Phylo.read(str(TREE), "newick")

tips = [t.name for t in tree.get_terminals() if t.name]

print("\nTree tips:", len(tips))

if len(tips) != 110:
    raise SystemExit(
        f"ERROR: expected 110 tree tips, found {len(tips)}"
    )

tip_to_tree_name = {
    extract_gcf(tip): tip
    for tip in tips
}

tree_gcf = set(
    x for x in tip_to_tree_name
    if x is not None
)

print("Tree GCF accessions:", len(tree_gcf))

# ============================================================
# BGC EVIDENCE
# ============================================================

df = pd.read_csv(EVIDENCE, sep="\t")

print("\nBGC evidence rows:", len(df))

if len(df) != 582:
    raise SystemExit(
        f"ERROR: expected 582 BGC rows, found {len(df)}"
    )

# ------------------------------------------------------------
# Authoritative c0.7 GCF identifier
#
# Your merged table contains:
#   Family_x
#   Family_y
#   GCF_family_x
#   GCF_family_y
#
# They represent the same family here.
#
# Use GCF_family_x because it is the c0.7 validated GCF
# field inherited from GCF_summary_c0.7_validated.tsv.
# ------------------------------------------------------------

if "GCF_family_x" in df.columns:
    gcf_col = "GCF_family_x"
elif "Family_x" in df.columns:
    gcf_col = "Family_x"
elif "GCF_family" in df.columns:
    gcf_col = "GCF_family"
elif "Family" in df.columns:
    gcf_col = "Family"
else:
    raise SystemExit(
        "ERROR: no GCF family column found."
    )

print("Using GCF column:", gcf_col)

df["GCF"] = df[gcf_col].astype(str)

df["genome_accession"] = df["accession"].astype(str)

print("Unique GCFs:", df["GCF"].nunique())
print("Unique genomes:", df["genome_accession"].nunique())

if df["GCF"].nunique() != 111:
    raise SystemExit(
        f"ERROR: expected 111 c0.7 GCFs, found {df['GCF'].nunique()}"
    )

# ============================================================
# ACCESSION VALIDATION
# ============================================================

bgc_genomes = set(df["genome_accession"])

missing = bgc_genomes - tree_gcf

print("\nBGC genome accessions:", len(bgc_genomes))
print("Missing from tree:", len(missing))

if missing:
    print("\n".join(sorted(missing)))
    raise SystemExit(
        "ERROR: some BGC genomes are absent from the final tree."
    )

print("PASS: all 110 BGC genomes map to the final tree.")

# ============================================================
# GCF → GENOMES
# ============================================================

gcf_to_genomes = (
    df[
        ["GCF", "genome_accession"]
    ]
    .drop_duplicates()
    .groupby("GCF")["genome_accession"]
    .apply(sorted)
    .to_dict()
)

# ============================================================
# PHYLOGENETIC ANALYSIS
# ============================================================

rows = []

for gcf, genomes in sorted(gcf_to_genomes.items()):

    genomes = sorted(set(genomes))

    tree_names = [
        tip_to_tree_name[g]
        for g in genomes
        if g in tip_to_tree_name
    ]

    n_genomes = len(genomes)

    if len(tree_names) == 0:
        continue

    if len(tree_names) == 1:

        mrca = tree.common_ancestor(tree_names[0])
        mrca_tip_count = len(mrca.get_terminals())

        mean_dist = 0.0
        median_dist = 0.0

    else:

        mrca = tree.common_ancestor(tree_names)
        mrca_tip_count = len(mrca.get_terminals())

        mean_dist = mean_pairwise_distance(
            tree,
            tree_names
        )

        median_dist = median_pairwise_distance(
            tree,
            tree_names
        )

    if mrca_tip_count > 0:
        mrca_enrichment = (
            n_genomes / mrca_tip_count
        )
    else:
        mrca_enrichment = np.nan

    prevalence = (
        n_genomes / 110 * 100
    )

    # --------------------------------------------------------
    # Descriptive operational classes
    # --------------------------------------------------------

    if n_genomes == 1:
        restriction = "single-genome"

    elif (
        n_genomes <= 10
        and mrca_enrichment >= 0.50
    ):
        restriction = "strongly lineage-restricted"

    elif (
        prevalence < 10
        and mrca_enrichment >= 0.25
    ):
        restriction = "lineage-restricted"

    elif mrca_enrichment >= 0.50:
        restriction = "phylogenetically concentrated"

    else:
        restriction = "phylogenetically dispersed"

    rows.append({
        "GCF": gcf,
        "tree_genome_count": n_genomes,
        "tree_prevalence_pct": prevalence,
        "MRCA_tip_count": mrca_tip_count,
        "MRCA_enrichment": mrca_enrichment,
        "mean_pairwise_distance": mean_dist,
        "median_pairwise_distance": median_dist,
        "phylogenetic_restriction_class": restriction,
        "tree_members": ";".join(genomes)
    })

phy = pd.DataFrame(rows)

# ============================================================
# ADD VALIDATED GCF SUMMARY
# ============================================================

summary_cols = [
    "GCF",
    "GCF_genomes",
    "GCF_prevalence_pct",
    "GCF_category",
    "GCF_products",
    "GCF_classes"
]

missing_summary = [
    c for c in summary_cols
    if c not in df.columns
]

if missing_summary:
    raise SystemExit(
        "ERROR: missing validated GCF summary columns: "
        + ", ".join(missing_summary)
    )

gcf_summary = (
    df[summary_cols]
    .drop_duplicates("GCF")
)

phy = phy.merge(
    gcf_summary,
    on="GCF",
    how="left",
    validate="one_to_one"
)

# ============================================================
# VALIDATION OF GCF PREVALENCE
# ============================================================

if not np.allclose(
    phy["tree_genome_count"],
    phy["GCF_genomes"],
    equal_nan=False
):
    print(
        "\nWARNING: tree-derived and validated GCF genome counts differ."
    )

    check = phy[
        phy["tree_genome_count"] != phy["GCF_genomes"]
    ]

    print(check.head(20).to_string(index=False))

else:
    print(
        "\nPASS: tree-derived GCF genome counts match "
        "validated GCF summary."
    )

# ============================================================
# SAVE GCF TABLE
# ============================================================

phy = phy.sort_values(
    [
        "MRCA_enrichment",
        "GCF_prevalence_pct"
    ],
    ascending=[
        False,
        True
    ]
).reset_index(drop=True)

phy.to_csv(
    OUT_GCF,
    sep="\t",
    index=False
)

# ============================================================
# MERGE TO BGC LEVEL
# ============================================================

bgc_phy_cols = [
    "GCF",
    "MRCA_tip_count",
    "MRCA_enrichment",
    "tree_fraction"
]

phy["tree_fraction"] = (
    phy["tree_genome_count"] / 110
)

bgc = df.merge(
    phy[
        bgc_phy_cols
        + [
            "mean_pairwise_distance",
            "median_pairwise_distance",
            "phylogenetic_restriction_class"
        ]
    ],
    on="GCF",
    how="left",
    validate="many_to_one"
)

bgc["phylogenetic_priority_flag"] = (
    bgc[
        "phylogenetic_restriction_class"
    ]
    .isin([
        "single-genome",
        "strongly lineage-restricted",
        "lineage-restricted"
    ])
    .astype(int)
)

bgc.to_csv(
    OUT_BGC,
    sep="\t",
    index=False
)

# ============================================================
# REPORT
# ============================================================

print("\n" + "=" * 100)
print("PHYLOGENETIC RESTRICTION SUMMARY")
print("=" * 100)

print("\nRestriction classes:")
print(
    phy[
        "phylogenetic_restriction_class"
    ].value_counts()
)

print("\nGCF prevalence:")
print(
    phy[
        "GCF_prevalence_pct"
    ].describe()
)

print("\nRestricted GCFs:")
restricted = phy[
    phy["phylogenetic_restriction_class"].isin([
        "single-genome",
        "strongly lineage-restricted",
        "lineage-restricted"
    ])
]

print(
    restricted[
        [
            "GCF",
            "GCF_genomes",
            "GCF_prevalence_pct",
            "MRCA_tip_count",
            "MRCA_enrichment",
            "mean_pairwise_distance",
            "phylogenetic_restriction_class",
            "GCF_classes",
            "GCF_products"
        ]
    ]
    .sort_values(
        [
            "MRCA_enrichment",
            "GCF_prevalence_pct"
        ],
        ascending=[
            False,
            True
        ]
    )
    .to_string(index=False)
)

print("\n" + "=" * 100)
print("FILES")
print("=" * 100)
print(OUT_GCF)
print(OUT_BGC)
