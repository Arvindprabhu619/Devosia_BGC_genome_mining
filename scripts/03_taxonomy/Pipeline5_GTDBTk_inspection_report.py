#!/usr/bin/env python3

# ============================================================
# Pipeline 5 — GTDB-Tk R232 classification inspection
#
# Purpose:
#   Inspect the GTDB-Tk bac120 classification summary,
#   extract taxonomic ranks, and summarize Devosia-related
#   and non-Devosia genomes.
#
# Input:
#   gtdbtk_r232/classify/gtdbtk.bac120.summary.tsv
#
# Output:
#   GTDBTk_R232_taxonomy_parsed.tsv
# ============================================================

from pathlib import Path

import pandas as pd


# ---- Determine repository root ----

PROJECT_ROOT = Path(__file__).resolve().parents[2]

GTDB_DIR = PROJECT_ROOT / "gtdbtk_r232"

INPUT_FILE = (
    GTDB_DIR
    / "classify"
    / "gtdbtk.bac120.summary.tsv"
)

OUTPUT_FILE = PROJECT_ROOT / "GTDBTk_R232_taxonomy_parsed.tsv"


# ---- Validate input ----

if not INPUT_FILE.exists():

    raise FileNotFoundError(
        f"GTDB-Tk summary file not found:\n"
        f"  {INPUT_FILE}\n\n"
        f"Run Pipeline4_GTDBTK_MAG.sh first."
    )


print("=" * 80)
print("GTDB-Tk R232 — COMPLETE BAC120 CLASSIFICATION INSPECTION")
print("=" * 80)

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)

print("\nTotal genomes classified:", len(df))
print("Total columns:", len(df.columns))


print("\n=== Classification methods ===")

print(
    df["classification_method"]
    .value_counts(dropna=False)
    .to_string()
)


print("\n=== Warnings ===")

print(
    df["warnings"]
    .value_counts(dropna=False)
    .to_string()
)


print("\n=== Translation tables ===")

print(
    df["translation_table"]
    .value_counts(dropna=False)
    .to_string()
)


# ------------------------------------------------------------
# Taxonomic extraction
# ------------------------------------------------------------

def extract_taxonomy(classification):

    ranks = {
        "domain": "d__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }

    result = {}

    if pd.isna(classification):

        for rank in ranks:
            result[rank] = "Unknown"

        return result

    parts = [
        x.strip()
        for x in str(classification).split(";")
    ]

    for rank, prefix in ranks.items():

        value = next(
            (
                x[len(prefix):]
                for x in parts
                if x.startswith(prefix)
            ),
            ""
        )

        result[rank] = (
            value
            if value
            else "Unknown"
        )

    return result


taxonomy = df["classification"].apply(
    extract_taxonomy
)

taxdf = pd.DataFrame(
    taxonomy.tolist()
)

for col in taxdf.columns:
    df[col] = taxdf[col]


print("\n=== GENUS COUNTS ===")

print(
    df["genus"]
    .value_counts()
    .to_string()
)


print("\n=== SPECIES COUNTS — TOP 30 ===")

print(
    df["species"]
    .value_counts()
    .head(30)
    .to_string()
)


print("\n=== DEVOSIA-RELATED GENOMES ===")

dev = df[
    df["genus"]
    .str.contains(
        "Devosia",
        case=False,
        na=False
    )
]

print(
    "Genomes with genus containing 'Devosia':",
    len(dev)
)

dev_columns = [
    "user_genome",
    "genus",
    "species",
    "classification_method",
    "closest_genome_reference",
    "closest_genome_ani",
    "closest_genome_af",
    "note",
]

print(
    dev[dev_columns]
    .to_string(index=False)
)


print("\n=== NON-DEVOSIA GENOMES ===")

nondev = df[
    ~df["genus"]
    .str.contains(
        "Devosia",
        case=False,
        na=False
    )
]

print(
    "Non-Devosia genomes:",
    len(nondev)
)

print(
    nondev[dev_columns]
    .to_string(index=False)
)


# ------------------------------------------------------------
# Save parsed taxonomy
# ------------------------------------------------------------

df.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUTPUT_FILE)
