
#!/usr/bin/env python3

import os
from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
BASE = SCRIPT_DIR.parent.parent.parent

MEMBER = os.path.join(
    BASE,
    "FAM_00321_member_BGCs_c0.7.tsv"
)

NOVELTY = os.path.join(
    BASE,
    "bgc_novelty_extended.tsv"
)

OUT = os.path.join(
    BASE,
    "FAM_00321_final_validation.tsv"
)

print("=" * 100)
print("FAM_00321 FINAL VALIDATION")
print("=" * 100)

# ------------------------------------------------------------
# Load files
# ------------------------------------------------------------

member = pd.read_csv(MEMBER, sep="\t", dtype=str)
novelty = pd.read_csv(NOVELTY, sep="\t", dtype=str)

print(f"\nFAM_00321 member BGCs : {len(member)}")
print(f"Novelty records       : {len(novelty)}")

# ------------------------------------------------------------
# Normalize identifiers
# ------------------------------------------------------------

for df in [member, novelty]:
    df["accession"] = df["accession"].astype(str).str.strip()
    df["contig"] = df["contig"].astype(str).str.strip()
    df["region_number"] = (
        pd.to_numeric(df["region_number"], errors="coerce")
        .astype("Int64")
    )

# ------------------------------------------------------------
# Check duplicate member identifiers
# ------------------------------------------------------------

member["member_key"] = (
    member["accession"] + "|" +
    member["contig"] + "|" +
    member["region_number"].astype(str)
)

novelty["member_key"] = (
    novelty["accession"] + "|" +
    novelty["contig"] + "|" +
    novelty["region_number"].astype(str)
)

print("\nDuplicate member keys:")
print(member["member_key"].duplicated().sum())

# ------------------------------------------------------------
# Merge novelty
# ------------------------------------------------------------

novelty_cols = [
    "member_key",
    "max_mibig_similarity",
    "top_mibig_hit",
    "top_mibig_description",
    "top_mibig_cluster_type",
    "top_hit_blast_score",
    "top_hit_n_homologous_genes",
    "top_hit_core_gene_hits",
    "novelty_category"
]

merged = member.merge(
    novelty[novelty_cols],
    on="member_key",
    how="left",
    validate="one_to_one"
)

# ------------------------------------------------------------
# Mapping QC
# ------------------------------------------------------------

mapped = merged["novelty_category"].notna().sum()
unmapped = merged["novelty_category"].isna().sum()

print("\n" + "=" * 100)
print("NOVELTY MAPPING QC")
print("=" * 100)

print(f"Total FAM_00321 members : {len(merged)}")
print(f"Mapped                  : {mapped}")
print(f"Unmapped                : {unmapped}")

if unmapped > 0:
    print("\nUNMAPPED MEMBERS:")
    print(
        merged.loc[
            merged["novelty_category"].isna(),
            ["GCF", "accession", "contig", "region_number", "GBK", "Record"]
        ].to_string(index=False)
    )

# ------------------------------------------------------------
# Biosynthetic identity
# ------------------------------------------------------------

print("\n" + "=" * 100)
print("BIOSYNTHETIC IDENTITY")
print("=" * 100)

print("\nProducts:")
print(merged["products"].value_counts(dropna=False).to_string())

print("\nBGC classes:")
print(merged["bgc_class"].value_counts(dropna=False).to_string())

# ------------------------------------------------------------
# Contig-edge distribution
# ------------------------------------------------------------

print("\n" + "=" * 100)
print("GENOMIC LOCATION")
print("=" * 100)

print("\nContig-edge status:")
print(merged["contig_edge"].value_counts(dropna=False).to_string())

# ------------------------------------------------------------
# Novelty distribution
# ------------------------------------------------------------

print("\n" + "=" * 100)
print("NOVELTY DISTRIBUTION")
print("=" * 100)

novelty_counts = (
    merged["novelty_category"]
    .value_counts(dropna=False)
    .rename_axis("novelty_category")
    .reset_index(name="n_BGCs")
)

novelty_counts["percentage"] = (
    novelty_counts["n_BGCs"] / len(merged) * 100
)

print(novelty_counts.to_string(index=False))

# ------------------------------------------------------------
# MIBiG similarity summary
# ------------------------------------------------------------

similarity = pd.to_numeric(
    merged["max_mibig_similarity"],
    errors="coerce"
)

print("\n" + "=" * 100)
print("MIBiG SIMILARITY")
print("=" * 100)

print(f"Non-zero MIBiG similarity : {(similarity > 0).sum()}")
print(f"Zero MIBiG similarity     : {(similarity == 0).sum()}")
print(f"Missing similarity        : {similarity.isna().sum()}")

if similarity.notna().any():
    print(f"Maximum similarity        : {similarity.max()}")
    print(f"Mean similarity           : {similarity.mean():.2f}")
    print(f"Median similarity         : {similarity.median():.2f}")

# ------------------------------------------------------------
# Number of genomes
# ------------------------------------------------------------

n_genomes = merged["accession"].nunique()

print("\n" + "=" * 100)
print("GENOME DISTRIBUTION")
print("=" * 100)

print(f"Unique genomes containing FAM_00321 : {n_genomes}")
print(f"Total analyzed genomes              : 110")
print(f"Prevalence                           : {n_genomes / 110 * 100:.2f}%")

# ------------------------------------------------------------
# Save final validation table
# ------------------------------------------------------------

merged.to_csv(
    OUT,
    sep="\t",
    index=False
)

print("\n" + "=" * 100)
print("OUTPUT")
print("=" * 100)

print(f"Written: {OUT}")

print("\nFINAL STATUS:")

if (
    len(merged) == 109
    and mapped == 109
    and n_genomes == 109
    and merged["products"].eq("terpene").all()
):
    print("PASS — FAM_00321 is completely validated.")
else:
    print("WARNING — FAM_00321 requires further inspection.")

print("=" * 100)

