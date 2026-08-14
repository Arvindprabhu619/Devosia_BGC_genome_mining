#!/usr/bin/env python3

import pandas as pd
import numpy as np

# ============================================================
# DEVOSIA BGC CANDIDATE PRIORITIZATION
#
# Phase 1:
#   Build 582-BGC master table
#   Integrate:
#       - BGC coordinates/products
#       - MIBiG similarity / novelty
#       - KCB similarity
#       - c0.7 GCF assignment
#       - GCF prevalence
#       - biosynthetic class
#
# Preliminary score:
#   Novelty              0-3
#   MIBiG relationship   0-2
#   GCF rarity           0-2
#   Product/class        0-1
#   ------------------------------------------------
#   Current maximum       8
#
# Domain architecture and phylogenetic restriction
# will be added in Phase 2.
# ============================================================

REGION = "bgc_regions_extended.tsv"
NOVELTY = "BGC_582_GTDB_R232_novelty_mapping.tsv"
KCB = "BGC_KnownClusterBlast_novelty_GTDB_R232.tsv"
GCF_MEM = "gcf_membership_full.tsv"
GCF_SUM = "GCF_summary_c0.7_validated.tsv"

OUT_MASTER = "BGC_candidate_master_c0.7.tsv"
OUT_TOP30 = "BGC_candidate_top30_c0.7.tsv"
OUT_SCORE_SUMMARY = "BGC_candidate_scoring_summary_c0.7.tsv"


def key(df):
    return (
        df["accession"].astype(str)
        + "|"
        + df["contig"].astype(str)
        + "|"
        + df["region_number"].astype(str)
    )


print("=" * 100)
print("DEVOSIA BGC CANDIDATE PRIORITIZATION — PHASE 1")
print("=" * 100)


# ------------------------------------------------------------
# 1. LOAD BGC REGIONS
# ------------------------------------------------------------

regions = pd.read_csv(REGION, sep="\t")

print("\nBGC regions:", regions.shape)

regions["BGC_key"] = key(regions)


# ------------------------------------------------------------
# 2. LOAD MIBiG / NOVELTY DATA
# ------------------------------------------------------------

nov = pd.read_csv(NOVELTY, sep="\t")

print("Novelty table:", nov.shape)

nov["BGC_key"] = key(nov)

nov_keep = [
    "BGC_key",
    "max_mibig_similarity",
    "top_mibig_hit",
    "top_mibig_description",
    "top_mibig_cluster_type",
    "top_hit_blast_score",
    "top_hit_n_homologous_genes",
    "top_hit_core_gene_hits",
    "novelty_category",
    "genus",
    "species",
]

nov = nov[nov_keep]


# ------------------------------------------------------------
# 3. LOAD KCB DATA
# ------------------------------------------------------------

kcb = pd.read_csv(KCB, sep="\t")

print("KCB table:", kcb.shape)

kcb["BGC_key"] = key(kcb)

kcb_keep = [
    "BGC_key",
    "has_kcb_result",
    "n_ranked_hits",
    "top_mibig",
    "top_mibig_type",
    "similarity_percent",
]

kcb = kcb[kcb_keep]


# ------------------------------------------------------------
# 4. LOAD GCF ASSIGNMENTS AT c0.7
# ------------------------------------------------------------

gcf = pd.read_csv(GCF_MEM, sep="\t")

gcf = gcf[gcf["cutoff"].astype(str) == "0.7"].copy()

print("c0.7 GCF membership rows:", len(gcf))

# Extract family and genomic record
#
# Record:
#   NZ_xxx.region001.gbk_region_1
#
# We reconstruct the BGC key from GBK / Record information.

gcf["Record"] = gcf["Record"].astype(str)

# Match region identifiers against the BGC-region structure.
#
# Extract accession from GBK:
# NZ_xxx.region001
#
# and region number from Record_Number.

gcf["accession"] = gcf["GBK"].astype(str).str.replace(
    r"\.region\d+$", "", regex=True
)

gcf["region_number"] = pd.to_numeric(
    gcf["Record_Number"], errors="coerce"
)

# GCF membership contains contig accession in GBK.
gcf["contig"] = gcf["accession"]

# We cannot use the above directly as BGC_key because the
# BGC-region table's accession corresponds to genome accession,
# while its contig is separate.
#
# Therefore create a region lookup using:
#   GBK = contig.regionXXX
#
# Extract contig from GBK.
gcf["contig"] = gcf["GBK"].astype(str).str.extract(
    r"^(.+)\.region\d+$"
)[0]

gcf["BGC_key_gcf"] = (
    gcf["contig"].astype(str)
    + "|"
    + gcf["region_number"].astype("Int64").astype(str)
)

# ------------------------------------------------------------
# 5. JOIN GCF MEMBERSHIP THROUGH CONTIG + REGION
# ------------------------------------------------------------

regions["region_lookup"] = (
    regions["contig"].astype(str)
    + "|"
    + regions["region_number"].astype(str)
)

gcf_lookup = gcf[
    [
        "region_lookup"
        if "region_lookup" in gcf.columns
        else "BGC_key_gcf",
        "Family",
        "GCF_ID",
        "bgc_class",
    ]
].copy()

# Correct lookup field
if "region_lookup" not in gcf.columns:
    gcf["region_lookup"] = gcf["contig"].astype(str) + "|" + \
        gcf["region_number"].astype("Int64").astype(str)

gcf_lookup = gcf[
    ["region_lookup", "Family", "GCF_ID", "bgc_class"]
].drop_duplicates("region_lookup")


# ------------------------------------------------------------
# 6. JOIN ALL DATA
# ------------------------------------------------------------

master = regions.merge(
    nov,
    on="BGC_key",
    how="left",
    validate="one_to_one"
)

master = master.merge(
    kcb,
    on="BGC_key",
    how="left",
    validate="one_to_one"
)

master = master.merge(
    gcf_lookup,
    on="region_lookup",
    how="left",
    validate="many_to_one"
)

print("\nAfter integration:")
print("Rows:", len(master))

if len(master) != 582:
    print("WARNING: expected 582 BGCs")


# ------------------------------------------------------------
# 7. JOIN GCF PREVALENCE
# ------------------------------------------------------------

gcf_sum = pd.read_csv(GCF_SUM, sep="\t")

gcf_sum = gcf_sum[
    [
        "GCF",
        "Genomes",
        "Prevalence_pct",
        "BGC_regions",
        "BGC_classes",
        "Products",
        "Category",
    ]
].copy()

gcf_sum = gcf_sum.rename(columns={
    "GCF": "GCF_family",
    "Genomes": "GCF_genomes",
    "Prevalence_pct": "GCF_prevalence_pct",
    "BGC_regions": "GCF_BGC_regions",
    "BGC_classes": "GCF_classes",
    "Products": "GCF_products",
    "Category": "GCF_category",
})

# Family is the clean GCF identifier
master["GCF_family"] = master["Family"]

master = master.merge(
    gcf_sum,
    on="GCF_family",
    how="left",
    validate="many_to_one"
)


# ------------------------------------------------------------
# 8. CLEAN MIBiG / NOVELTY VALUES
# ------------------------------------------------------------

master["max_mibig_similarity"] = pd.to_numeric(
    master["max_mibig_similarity"],
    errors="coerce"
).fillna(0)

master["similarity_percent"] = pd.to_numeric(
    master["similarity_percent"],
    errors="coerce"
)

master["similarity_percent"] = master["similarity_percent"].fillna(
    master["max_mibig_similarity"]
)

master["novelty_category"] = master["novelty_category"].fillna(
    "unknown"
)


# ------------------------------------------------------------
# 9. NOVELTY SCORE
# ------------------------------------------------------------

def novelty_score(x):
    if x == "novel" or x == "putatively_novel":
        return 3
    elif x == "related":
        return 1
    elif x == "known":
        return 0
    return 0


master["score_novelty"] = master["novelty_category"].apply(
    novelty_score
)


# ------------------------------------------------------------
# 10. MIBiG SCORE
#
# We don't simply reward high similarity.
#
# The scoring emphasizes either:
#
#   A) strong novelty (very low similarity)
#   B) potentially informative intermediate similarity
#
# while known/highly similar clusters receive lower
# discovery-priority scores.
# ------------------------------------------------------------

def mibig_score(sim):
    if pd.isna(sim) or sim == 0:
        return 2
    elif sim <= 10:
        return 2
    elif sim <= 40:
        return 2
    elif sim <= 70:
        return 1
    else:
        return 0


master["score_mibig"] = master["similarity_percent"].apply(
    mibig_score
)


# ------------------------------------------------------------
# 11. GCF RARITY SCORE
# ------------------------------------------------------------

def rarity_score(p):
    if pd.isna(p):
        return 0
    elif p < 10:
        return 2
    elif p < 90:
        return 1
    else:
        return 0


master["score_gcf_rarity"] = master["GCF_prevalence_pct"].apply(
    rarity_score
)


# ------------------------------------------------------------
# 12. BIOSYNTHETIC CLASS / PRODUCT SCORE
#
# High-priority discovery classes:
#   PKS
#   NRPS
#   PKS/NRPS hybrids
#   RiPP
#   phosphonate
#   metallophore-associated systems
#
# We deliberately give this only one point so that
# class alone cannot dominate the ranking.
# ------------------------------------------------------------

def class_score(row):

    text = " ".join([
        str(row.get("products", "")),
        str(row.get("bgc_class", "")),
        str(row.get("GCF_classes", "")),
    ]).lower()

    high_interest = [
        "pks",
        "nrps",
        "ripp",
        "phosphonate",
        "metallophore",
        "lanthipeptide",
        "lassopeptide",
        "thioamitide",
    ]

    if any(x in text for x in high_interest):
        return 1

    return 0


master["score_product_class"] = master.apply(
    class_score,
    axis=1
)


# ------------------------------------------------------------
# 13. PRELIMINARY TOTAL
# ------------------------------------------------------------

master["score_preliminary"] = (
    master["score_novelty"]
    + master["score_mibig"]
    + master["score_gcf_rarity"]
    + master["score_product_class"]
)

# Current maximum = 8


# ------------------------------------------------------------
# 14. PRIORITY CATEGORY
# ------------------------------------------------------------

def priority_category(score):

    if score >= 7:
        return "High"
    elif score >= 5:
        return "Medium-high"
    elif score >= 3:
        return "Moderate"
    else:
        return "Lower"


master["priority_category"] = master[
    "score_preliminary"
].apply(priority_category)


# ------------------------------------------------------------
# 15. BIOLOGICAL RATIONALE
# ------------------------------------------------------------

def rationale(row):

    reasons = []

    if row["score_novelty"] >= 3:
        reasons.append("putatively novel BGC")

    if row["score_mibig"] >= 2:
        reasons.append("low/no MIBiG similarity")

    if row["score_gcf_rarity"] == 2:
        reasons.append("rare GCF")

    elif row["score_gcf_rarity"] == 1:
        reasons.append("intermediate-prevalence GCF")

    if row["score_product_class"] == 1:
        reasons.append("biosynthetically interesting product class")

    if not reasons:
        reasons.append("comparative biosynthetic interest")

    return "; ".join(reasons)


master["priority_rationale"] = master.apply(
    rationale,
    axis=1
)


# ------------------------------------------------------------
# 16. SORT
# ------------------------------------------------------------

master = master.sort_values(
    [
        "score_preliminary",
        "score_novelty",
        "score_gcf_rarity",
        "score_mibig",
    ],
    ascending=[False, False, False, False]
).reset_index(drop=True)

master["preliminary_rank"] = np.arange(1, len(master) + 1)


# ------------------------------------------------------------
# 17. SAVE MASTER
# ------------------------------------------------------------

master.to_csv(
    OUT_MASTER,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# 18. TOP 30
# ------------------------------------------------------------

top30 = master.head(30).copy()

top30.to_csv(
    OUT_TOP30,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# 19. SCORING SUMMARY
# ------------------------------------------------------------

summary = pd.DataFrame({
    "criterion": [
        "Novelty",
        "MIBiG relationship",
        "GCF rarity",
        "Product/class",
        "Domain architecture",
        "Phylogenetic restriction",
    ],
    "maximum_points_current_phase": [
        3,
        2,
        2,
        1,
        0,
        0,
    ],
    "status": [
        "included",
        "included",
        "included",
        "included",
        "Phase 2",
        "Phase 2",
    ],
    "rationale": [
        "Prioritizes putatively novel biosynthetic potential",
        "Distinguishes unexplored/divergent clusters from strongly characterized clusters",
        "Prioritizes lineage-restricted GCFs",
        "Modestly prioritizes PKS/NRPS/RiPP and unusual biosynthetic classes",
        "Requires antiSMASH domain/architecture inspection",
        "Requires integration with final Devosia phylogeny",
    ]
})

summary.to_csv(
    OUT_SCORE_SUMMARY,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# 20. REPORT
# ------------------------------------------------------------

print("\n" + "=" * 100)
print("VALIDATION")
print("=" * 100)

print("\nMaster rows:", len(master))
print("Unique BGC keys:", master["BGC_key"].nunique())

print("\nNovelty:")
print(master["novelty_category"].value_counts())

print("\nGCF assignment:")
print("Assigned:", master["GCF_family"].notna().sum())
print("Missing:", master["GCF_family"].isna().sum())

print("\nGCF categories:")
print(master["GCF_category"].value_counts(dropna=False))

print("\nPreliminary score:")
print(
    master["score_preliminary"]
    .value_counts()
    .sort_index(ascending=False)
)

print("\nPriority category:")
print(master["priority_category"].value_counts())

print("\n" + "=" * 100)
print("TOP 30 PRELIMINARY CANDIDATES")
print("=" * 100)

cols = [
    "preliminary_rank",
    "accession",
    "contig",
    "region_number",
    "start",
    "end",
    "products",
    "novelty_category",
    "max_mibig_similarity",
    "similarity_percent",
    "Family",
    "GCF_genomes",
    "GCF_prevalence_pct",
    "GCF_category",
    "score_novelty",
    "score_mibig",
    "score_gcf_rarity",
    "score_product_class",
    "score_preliminary",
    "priority_category",
    "priority_rationale",
]

print(
    top30[cols].to_string(index=False)
)

print("\n" + "=" * 100)
print("FILES GENERATED")
print("=" * 100)

print(OUT_MASTER)
print(OUT_TOP30)
print(OUT_SCORE_SUMMARY)
