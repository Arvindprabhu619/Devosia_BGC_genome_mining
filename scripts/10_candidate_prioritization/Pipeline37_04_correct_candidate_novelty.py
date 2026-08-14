#!/usr/bin/env python3

import pandas as pd

MASTER = "BGC_candidate_master_c0.7.tsv"
KCB = "BGC_KnownClusterBlast_novelty_GTDB_R232.tsv"

OUT = "BGC_candidate_master_c0.7_corrected.tsv"

master = pd.read_csv(
    MASTER,
    sep="\t"
)

kcb = pd.read_csv(
    KCB,
    sep="\t"
)

keys = [
    "accession",
    "contig",
    "region_number"
]

kcb = kcb[
    keys +
    [
        "has_kcb_result",
        "n_ranked_hits",
        "top_mibig",
        "top_mibig_type",
        "similarity_percent",
        "novelty_category"
    ]
].copy()

kcb = kcb.rename(columns={
    "similarity_percent": "KCB_similarity_percent",
    "novelty_category": "KCB_novelty_category"
})

out = master.drop(
    columns=[
        "novelty_category",
        "similarity_percent"
    ],
    errors="ignore"
).merge(
    kcb,
    on=keys,
    how="left",
    validate="one_to_one"
)

# ------------------------------------------------------------
# Authoritative novelty classification
# Matches manuscript methodology:
#
# >=70% = known
# >0 and <70% = related
# no KCB result = putatively novel
# ------------------------------------------------------------

def classify(row):

    sim = row["KCB_similarity_percent"]

    if pd.isna(sim):
        return "putatively_novel"

    if sim >= 70:
        return "known"

    return "related"


out["novelty_category"] = out.apply(
    classify,
    axis=1
)

out["similarity_percent"] = out[
    "KCB_similarity_percent"
]

# ------------------------------------------------------------
# Recalculate novelty score
# ------------------------------------------------------------

out["score_novelty"] = (
    out["novelty_category"]
    .map({
        "putatively_novel": 3,
        "related": 1,
        "known": 0
    })
    .fillna(0)
    .astype(int)
)

# ------------------------------------------------------------
# Recalculate MIBiG score using KCB similarity
#
# 0%/no match -> novelty evidence
# 1-40%       -> divergent related space
# 41-70%      -> moderately related
# >70%        -> known
# ------------------------------------------------------------

def mibig_score(sim):

    if pd.isna(sim) or sim == 0:
        return 2
    elif sim <= 40:
        return 2
    elif sim < 70:
        return 1
    else:
        return 0


out["score_mibig"] = (
    out["similarity_percent"]
    .apply(mibig_score)
)

out["score_preliminary"] = (
    out["score_novelty"]
    + out["score_mibig"]
    + out["score_gcf_rarity"]
    + out["score_product_class"]
)

out.to_csv(
    OUT,
    sep="\t",
    index=False
)

print("="*100)
print("CORRECTED NOVELTY CLASSIFICATION")
print("="*100)

print(
    out["novelty_category"]
    .value_counts()
)

print("\nSimilarity:")
print(
    out["similarity_percent"].describe()
)

print("\nExpected manuscript classification:")
print("Novel  :", (out["novelty_category"] == "putatively_novel").sum())
print("Related:", (out["novelty_category"] == "related").sum())
print("Known  :", (out["novelty_category"] == "known").sum())

assert (
    out["novelty_category"].value_counts().to_dict()
    == {
        "putatively_novel": 476,
        "related": 70,
        "known": 36
    }
), "Novelty counts do not match final manuscript classification."

print("\nPASS: 476 / 70 / 36")

print("\nSaved:")
print(OUT)
