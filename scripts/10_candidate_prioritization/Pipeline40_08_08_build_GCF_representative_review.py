#!/usr/bin/env python3

import pandas as pd
import numpy as np

INPUT = "BGC_integrated_candidate_pool_c0.7.tsv"
KCB = "BGC_KnownClusterBlast_novelty_GTDB_R232.tsv"

OUT = "GCF_representative_candidate_review_c0.7.tsv"

df = pd.read_csv(INPUT, sep="\t")
kcb = pd.read_csv(KCB, sep="\t")

KEYS = [
    "accession",
    "contig",
    "region_number"
]

# ------------------------------------------------------------
# Add authoritative KCB information
# ------------------------------------------------------------

kcb_keep = kcb[
    KEYS + [
        "top_mibig",
        "top_mibig_type",
        "similarity_percent",
        "n_ranked_hits",
        "note"
    ]
].copy()

kcb_keep = kcb_keep.rename(columns={
    "top_mibig": "MIBiG_hit",
    "top_mibig_type": "MIBiG_type",
    "similarity_percent": "MIBiG_similarity",
    "n_ranked_hits": "MIBiG_ranked_hits",
    "note": "KCB_note"
})

df = df.drop(
    columns=[
        "MIBiG_hit",
        "MIBiG_similarity",
        "MIBiG_ranked_hits"
    ],
    errors="ignore"
).merge(
    kcb_keep,
    on=KEYS,
    how="left",
    validate="one_to_one"
)

# ------------------------------------------------------------
# Normalize
# ------------------------------------------------------------

for c in [
    "MIBiG_similarity",
    "GCF_prevalence_final",
    "candidate_evidence_points",
    "architecture_score_final",
    "aSDomain_count_final",
    "PFAM_count_final",
    "GCF_genomes_final",
    "MRCA_enrichment_final",
]:
    if c in df.columns:
        df[c] = pd.to_numeric(
            df[c],
            errors="coerce"
        )

df["MIBiG_similarity"] = df[
    "MIBiG_similarity"
].fillna(0)

df["edge_flag"] = (
    df["assembly_priority"]
    == "edge-region-review"
)

# ------------------------------------------------------------
# Representative BGC selection within each GCF
#
# Preference:
#   1. non-edge
#   2. highest evidence
#   3. distinctive architecture
#   4. higher domain evidence
#   5. lower MIBiG similarity for novel clusters
# ------------------------------------------------------------

df["representative_rank"] = 0.0

df["representative_rank"] += (
    (~df["edge_flag"]).astype(int) * 100
)

df["representative_rank"] += (
    df["candidate_evidence_points"]
    .fillna(0) * 10
)

df["representative_rank"] += (
    df["architecture_score_final"]
    .fillna(0) * 5
)

df["representative_rank"] += (
    df["aSDomain_count_final"]
    .fillna(0) / 100
)

# For novel BGCs, prefer zero/no MIBiG similarity
df.loc[
    df["novelty_final"] == "putatively_novel",
    "representative_rank"
] += (
    (100 - df["MIBiG_similarity"])
    / 100
)

# ------------------------------------------------------------
# One representative per GCF
# ------------------------------------------------------------

df = df.sort_values(
    [
        "GCF_final",
        "representative_rank"
    ],
    ascending=[
        True,
        False
    ]
)

representatives = (
    df.groupby(
        "GCF_final",
        as_index=False,
        group_keys=False
    )
    .head(1)
    .copy()
)

print("="*110)
print("GCF REPRESENTATIVE CANDIDATE REVIEW")
print("="*110)

print(
    "Candidate BGCs:",
    len(df)
)

print(
    "Candidate GCFs:",
    df["GCF_final"].nunique()
)

print(
    "One representative per GCF:",
    len(representatives)
)

# ------------------------------------------------------------
# Biological candidate class
# ------------------------------------------------------------

def candidate_type(row):

    if (
        row["novelty_final"] == "putatively_novel"
        and
        row["rare_gcf"]
        and
        row["distinctive_architecture"]
        and
        row["multi_genome_lineage_restricted"]
    ):
        return "Novel + rare + unusual + lineage-restricted"

    if (
        row["novelty_final"] == "putatively_novel"
        and
        row["rare_gcf"]
        and
        row["distinctive_architecture"]
    ):
        return "Novel + rare + unusual architecture"

    if (
        row["novelty_final"] == "putatively_novel"
        and
        row["multi_genome_lineage_restricted"]
    ):
        return "Novel + lineage-restricted"

    if (
        row["novelty_final"] == "related"
        and
        row["distinctive_architecture"]
    ):
        return "Divergent known-space + unusual architecture"

    if row["single_genome"]:
        return "Singleton / highly restricted"

    return "Other"

representatives["candidate_type"] = representatives.apply(
    candidate_type,
    axis=1
)

# ------------------------------------------------------------
# Add a chemistry grouping
# ------------------------------------------------------------

def chemistry_group(p):

    p = str(p).lower()

    if "pks" in p and "nrps" in p:
        return "PKS-NRPS"

    if "pks" in p:
        return "PKS"

    if any(
        x in p
        for x in [
            "ripp",
            "lassopeptide",
            "thioamitide",
            "thiopeptide",
            "lanthipeptide",
            "lap",
            "rre-containing"
        ]
    ):
        return "RiPP"

    if "nrps" in p or "metallophore" in p:
        return "NRPS/metallophore"

    if "phosphonate" in p:
        return "phosphonate"

    if "terpene" in p:
        return "terpene"

    if "betalactone" in p:
        return "betalactone"

    return "other"

representatives["chemistry_group"] = representatives[
    "product_final"
].apply(chemistry_group)

# ------------------------------------------------------------
# Compact review table
# ------------------------------------------------------------

review_cols = [
    "GCF_final",
    "accession",
    "contig",
    "region_number",
    "start",
    "end",
    "product_final",
    "chemistry_group",
    "novelty_final",
    "MIBiG_similarity",
    "MIBiG_hit",
    "MIBiG_type",
    "GCF_genomes_final",
    "GCF_prevalence_final",
    "GCF_category",
    "architecture_category",
    "aSDomain_count_final",
    "PFAM_count_final",
    "n_candidate_clusters_x"
]

review_cols = [
    c for c in review_cols
    if c in representatives.columns
]

review = representatives[
    review_cols
    + [
        "MRCA_enrichment_final",
        "phylogenetic_restriction_class",
        "candidate_evidence_points",
        "candidate_type",
        "assembly_priority"
    ]
].copy()

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

review.to_csv(
    OUT,
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# Report grouped by chemistry
# ------------------------------------------------------------

print("\n" + "="*110)
print("REPRESENTATIVE GCFs BY CHEMISTRY")
print("="*110)

print(
    review[
        "chemistry_group"
    ].value_counts()
)

print("\n" + "="*110)
print("REPRESENTATIVE GCF CANDIDATES")
print("="*110)

print(
    review.sort_values(
        [
            "candidate_evidence_points",
            "GCF_prevalence_final"
        ],
        ascending=[
            False,
            True
        ]
    ).to_string(
        index=False
    )
)

print("\nSaved:")
print(OUT)
