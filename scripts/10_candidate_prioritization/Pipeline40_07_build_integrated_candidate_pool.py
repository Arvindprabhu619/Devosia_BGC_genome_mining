#!/usr/bin/env python3

import pandas as pd
import numpy as np

INPUT = "BGC_candidate_evidence_phylogeny_c0.7.tsv"
OUT_POOL = "BGC_integrated_candidate_pool_c0.7.tsv"


def first_existing(df, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c

    if required:
        print("\nAVAILABLE COLUMNS:")
        print("\n".join(df.columns))
        raise SystemExit(
            "\nERROR: none of these columns found:\n  "
            + "\n  ".join(candidates)
        )

    return None


print("=" * 110)
print("INTEGRATED DEVOSIA BGC CANDIDATE POOL")
print("=" * 110)

df = pd.read_csv(INPUT, sep="\t")

print("Input BGCs:", len(df))

if len(df) != 582:
    raise SystemExit(
        f"ERROR: expected 582 BGCs, found {len(df)}"
    )

# ============================================================
# Resolve columns explicitly
# ============================================================

PRODUCT_COL = first_existing(
    df,
    [
        "products_x",
        "products",
        "products_y",
    ]
)

GCF_COL = first_existing(
    df,
    [
        "GCF",
        "GCF_family_x",
        "Family_x",
        "GCF_family_y",
        "Family_y",
    ]
)

PREV_COL = first_existing(
    df,
    [
        "GCF_prevalence_pct",
        "GCF_prevalence_pct_x",
        "GCF_prevalence_pct_y",
    ]
)

GENOME_COL = first_existing(
    df,
    [
        "GCF_genomes",
        "GCF_genomes_x",
        "GCF_genomes_y",
    ]
)

CATEGORY_COL = first_existing(
    df,
    [
        "GCF_category",
        "GCF_category_x",
        "GCF_category_y",
    ]
)

NOVELTY_COL = first_existing(
    df,
    [
        "novelty_category",
        "KCB_novelty_category",
    ]
)

ARCH_COL = first_existing(
    df,
    [
        "architecture_category"
    ]
)

ARCH_SCORE_COL = first_existing(
    df,
    [
        "architecture_evidence_score"
    ]
)

MRCA_ENRICH_COL = first_existing(
    df,
    [
        "MRCA_enrichment"
    ]
)

MRCA_TIPS_COL = first_existing(
    df,
    [
        "MRCA_tip_count"
    ]
)

PHYLO_COL = first_existing(
    df,
    [
        "phylogenetic_restriction_class"
    ]
)

SIM_COL = first_existing(
    df,
    [
        "similarity_percent",
        "KCB_similarity_percent",
        "max_mibig_similarity",
    ]
)

CC_COL = first_existing(
    df,
    [
        "n_candidate_clusters_x",
        "n_candidate_clusters",
        "n_candidate_clusters_y",
    ]
)

ASDOMAIN_COL = first_existing(
    df,
    [
        "aSDomain_count"
    ]
)

PFAM_COL = first_existing(
    df,
    [
        "PFAM_feature_count"
    ]
)

EDGE_COL = first_existing(
    df,
    [
        "contig_edge_x",
        "contig_edge",
        "contig_edge_y",
    ],
    required=False
)

print("\nRESOLVED FIELDS")
print("=" * 80)

for label, col in [
    ("Physical BGC product", PRODUCT_COL),
    ("GCF identifier", GCF_COL),
    ("GCF prevalence", PREV_COL),
    ("GCF genome count", GENOME_COL),
    ("GCF category", CATEGORY_COL),
    ("Novelty", NOVELTY_COL),
    ("Architecture", ARCH_COL),
    ("Architecture evidence", ARCH_SCORE_COL),
    ("MRCA enrichment", MRCA_ENRICH_COL),
    ("MRCA tips", MRCA_TIPS_COL),
    ("Phylogenetic class", PHYLO_COL),
    ("MIBiG similarity", SIM_COL),
    ("Candidate-cluster count", CC_COL),
    ("aSDomain count", ASDOMAIN_COL),
    ("PFAM count", PFAM_COL),
    ("Contig edge", EDGE_COL),
]:
    print(f"{label:30s} -> {col}")

# ============================================================
# STANDARDIZED WORKING VARIABLES
# ============================================================

df["product_final"] = (
    df[PRODUCT_COL]
    .fillna("")
    .astype(str)
)

df["GCF_final"] = (
    df[GCF_COL]
    .fillna("")
    .astype(str)
)

df["novelty_final"] = (
    df[NOVELTY_COL]
    .fillna("")
    .astype(str)
)

df["GCF_prevalence_final"] = pd.to_numeric(
    df[PREV_COL],
    errors="coerce"
)

df["GCF_genomes_final"] = pd.to_numeric(
    df[GENOME_COL],
    errors="coerce"
)

df["architecture_score_final"] = pd.to_numeric(
    df[ARCH_SCORE_COL],
    errors="coerce"
).fillna(0)

df["MRCA_enrichment_final"] = pd.to_numeric(
    df[MRCA_ENRICH_COL],
    errors="coerce"
)

df["MRCA_tip_count_final"] = pd.to_numeric(
    df[MRCA_TIPS_COL],
    errors="coerce"
)

df["MIBiG_similarity_final"] = pd.to_numeric(
    df[SIM_COL],
    errors="coerce"
).fillna(0)

df["candidate_clusters_final"] = pd.to_numeric(
    df[CC_COL],
    errors="coerce"
).fillna(1)

df["aSDomain_count_final"] = pd.to_numeric(
    df[ASDOMAIN_COL],
    errors="coerce"
).fillna(0)

df["PFAM_count_final"] = pd.to_numeric(
    df[PFAM_COL],
    errors="coerce"
).fillna(0)

if EDGE_COL is not None:
    df["contig_edge_final"] = (
        df[EDGE_COL]
        .astype(str)
        .str.lower()
        .isin(["true", "1"])
    )
else:
    df["contig_edge_final"] = False

# ============================================================
# VALIDATION
# ============================================================

print("\n" + "=" * 110)
print("VALIDATION")
print("=" * 110)

print("\nNovelty:")
print(df["novelty_final"].value_counts())

expected_novelty = {
    "putatively_novel": 476,
    "related": 70,
    "known": 36
}

observed = df["novelty_final"].value_counts().to_dict()

if observed != expected_novelty:
    raise SystemExit(
        f"ERROR: novelty mismatch.\n"
        f"Expected: {expected_novelty}\n"
        f"Observed: {observed}"
    )

print("PASS: 476 / 70 / 36")

print(
    "\nUnique c0.7 GCFs:",
    df["GCF_final"].nunique()
)

if df["GCF_final"].nunique() != 111:
    raise SystemExit(
        "ERROR: expected 111 c0.7 GCFs."
    )

print("PASS: 111 GCFs")

# ============================================================
# BIOLOGICAL FLAGS
# ============================================================

df["novel_flag"] = (
    df["novelty_final"]
    == "putatively_novel"
)

df["related_flag"] = (
    df["novelty_final"]
    == "related"
)

df["rare_gcf"] = (
    df["GCF_prevalence_final"] < 10
)

df["intermediate_gcf"] = (
    (df["GCF_prevalence_final"] >= 10)
    &
    (df["GCF_prevalence_final"] < 90)
)

df["core_gcf"] = (
    df["GCF_prevalence_final"] >= 90
)

# ============================================================
# PRODUCT / CLASS INTEREST
# ============================================================

df["high_interest_class"] = (
    df["product_final"]
    .str.lower()
    .str.contains(
        "nrps|pks|ripp|"
        "thioamitide|thiopeptide|"
        "lanthipeptide|lassopeptide|"
        "phosphonate|metallophore",
        regex=True
    )
)

df["hybrid_product"] = (
    df["product_final"]
    .str.contains(";", regex=False)
)

# ============================================================
# DISTINCTIVE ARCHITECTURE
# ============================================================

df["distinctive_architecture"] = (
    df[ARCH_COL]
    .fillna("")
    .astype(str)
    .isin([
        "PKS-NRPS hybrid",
        "RRE/RiPP-associated",
        "multi-product/hybrid",
        "high-domain-complexity",
    ])
)

# ============================================================
# PHYLOGENETIC RESTRICTION
# ============================================================

df["multi_genome_lineage_restricted"] = (
    df[PHYLO_COL]
    .fillna("")
    .isin([
        "strongly lineage-restricted",
        "lineage-restricted",
    ])
    &
    (df["GCF_genomes_final"] >= 2)
)

df["single_genome"] = (
    df["GCF_genomes_final"] == 1
)

# ============================================================
# EVIDENCE SCREEN
# ============================================================

df["candidate_evidence_points"] = 0

# Novelty
df.loc[
    df["novel_flag"],
    "candidate_evidence_points"
] += 3

df.loc[
    df["related_flag"],
    "candidate_evidence_points"
] += 1

# Rare family
df.loc[
    df["rare_gcf"],
    "candidate_evidence_points"
] += 2

# Discovery-relevant class
df.loc[
    df["high_interest_class"],
    "candidate_evidence_points"
] += 1

# Hybrid product
df.loc[
    df["hybrid_product"],
    "candidate_evidence_points"
] += 1

# Distinctive architecture
df.loc[
    df["distinctive_architecture"],
    "candidate_evidence_points"
] += 2

# True multi-genome phylogenetic restriction
df.loc[
    df["multi_genome_lineage_restricted"],
    "candidate_evidence_points"
] += 2

# Very high aSDomain complexity
df.loc[
    df["aSDomain_count_final"] >= 18,
    "candidate_evidence_points"
] += 1

# Multiple candidate clusters
df.loc[
    df["candidate_clusters_final"] >= 2,
    "candidate_evidence_points"
] += 1

# ============================================================
# EVIDENCE CLASS
# ============================================================

def evidence_class(row):

    if (
        row["novel_flag"]
        and row["rare_gcf"]
        and row["distinctive_architecture"]
    ):
        return (
            "Novel + rare + distinctive architecture"
        )

    if (
        row["novel_flag"]
        and row["multi_genome_lineage_restricted"]
    ):
        return (
            "Novel + lineage-restricted"
        )

    if (
        row["related_flag"]
        and row["rare_gcf"]
        and row["distinctive_architecture"]
    ):
        return (
            "Divergent known-space + "
            "rare + distinctive architecture"
        )

    if (
        row["novel_flag"]
        and row["single_genome"]
        and row["high_interest_class"]
    ):
        return (
            "Rare singleton + "
            "discovery-relevant class"
        )

    if (
        row["hybrid_product"]
        and row["rare_gcf"]
    ):
        return "Rare hybrid-product BGC"

    return "Other"

df["candidate_evidence_class"] = df.apply(
    evidence_class,
    axis=1
)

# ============================================================
# CANDIDATE POOL
# ============================================================

keep = (
    (df["candidate_evidence_points"] >= 6)
    |
    (
        df["multi_genome_lineage_restricted"]
        &
        df["novel_flag"]
    )
    |
    (
        df["novel_flag"]
        &
        df["high_interest_class"]
        &
        df["rare_gcf"]
    )
    |
    (
        df["novel_flag"]
        &
        df["hybrid_product"]
        &
        df["rare_gcf"]
    )
)

pool = df.loc[keep].copy()

pool["assembly_priority"] = np.where(
    df.loc[pool.index, "contig_edge_final"],
    "edge-region-review",
    "preferred"
)

# ============================================================
# SORT
# ============================================================

pool = pool.sort_values(
    [
        "candidate_evidence_points",
        "novel_flag",
        "multi_genome_lineage_restricted",
        "rare_gcf",
        "distinctive_architecture",
        "architecture_score_final",
        "GCF_prevalence_final",
    ],
    ascending=[
        False,
        False,
        False,
        False,
        False,
        False,
        True
    ]
).reset_index(drop=True)

pool["pool_rank"] = np.arange(
    1,
    len(pool) + 1
)

# ============================================================
# SAVE
# ============================================================

pool.to_csv(
    OUT_POOL,
    sep="\t",
    index=False
)

# ============================================================
# REPORT
# ============================================================

print("\n" + "=" * 110)
print("CANDIDATE POOL SUMMARY")
print("=" * 110)

print("Pool size:", len(pool))

print(
    "Unique GCFs:",
    pool["GCF_final"].nunique()
)

print("\nEvidence classes:")
print(
    pool["candidate_evidence_class"]
    .value_counts()
)

print("\nNovelty:")
print(
    pool["novelty_final"]
    .value_counts()
)

print("\nArchitecture:")
print(
    pool[ARCH_COL]
    .value_counts()
)

print("\nPhylogenetic restriction:")
print(
    pool[PHYLO_COL]
    .value_counts()
)

print("\n" + "=" * 110)
print("TOP 60 BGC CANDIDATES")
print("=" * 110)

display_cols = [
    "pool_rank",
    "accession",
    "contig",
    "region_number",
    "product_final",
    "MIBiG_similarity_final",
    "novelty_final",
    "GCF_final",
    "GCF_genomes_final",
    "GCF_prevalence_final",
    "CATEGORY_COL_PLACEHOLDER",
    "ARCH_COL_PLACEHOLDER",
    "PHYLO_COL_PLACEHOLDER",
    "MRCA_enrichment_final",
    "architecture_score_final",
    "candidate_evidence_points",
    "candidate_evidence_class",
    "assembly_priority",
]

# Replace placeholders with actual columns
display_cols = [
    (
        CATEGORY_COL
        if c == "CATEGORY_COL_PLACEHOLDER"
        else c
    )
    for c in display_cols
]

display_cols = [
    (
        ARCH_COL
        if c == "ARCH_COL_PLACEHOLDER"
        else c
    )
    for c in display_cols
]

display_cols = [
    (
        PHYLO_COL
        if c == "PHYLO_COL_PLACEHOLDER"
        else c
    )
    for c in display_cols
]

print(
    pool[
        display_cols
    ]
    .head(60)
    .to_string(index=False)
)

# ============================================================
# GCF SUMMARY
# ============================================================

print("\n" + "=" * 110)
print("TOP CANDIDATE GCFs")
print("=" * 110)

gcf_summary = (
    pool.groupby("GCF_final")
    .agg(
        candidate_BGCs=(
            "accession",
            "count"
        ),
        novelty=(
            "novelty_final",
            lambda x: ";".join(
                sorted(set(x))
            )
        ),
        prevalence=(
            "GCF_prevalence_final",
            "first"
        ),
        GCF_category=(
            CATEGORY_COL,
            "first"
        ),
        products=(
            "product_final",
            lambda x: ";".join(
                sorted(set(x))
            )
        ),
        architecture=(
            ARCH_COL,
            lambda x: ";".join(
                sorted(set(x))
            )
        ),
        phylogeny=(
            PHYLO_COL,
            lambda x: ";".join(
                sorted(set(x))
            )
        ),
        max_evidence=(
            "candidate_evidence_points",
            "max"
        ),
    )
    .reset_index()
    .sort_values(
        [
            "max_evidence",
            "prevalence"
        ],
        ascending=[
            False,
            True
        ]
    )
)

print(
    gcf_summary
    .head(60)
    .to_string(index=False)
)

print("\nSaved:")
print(OUT_POOL)
