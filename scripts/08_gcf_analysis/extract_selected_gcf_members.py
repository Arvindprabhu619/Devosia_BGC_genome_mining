#!/usr/bin/env python3

import os
import glob
import re
import pandas as pd

BASE = "bigscape_output_multicutoff/output_files"
CUTOFF = "0.7"

SELECTED = {
    "FAM_00321",
    "FAM_00355",
    "FAM_00260",
}

PATTERN = os.path.join(
    BASE,
    f"*_c{CUTOFF}",
    "*",
    f"*_clustering_c{CUTOFF}.tsv"
)

files = sorted(glob.glob(PATTERN))

print(f"Found {len(files)} clustering files")

frames = []

for path in files:

    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        print(f"[WARNING] Could not read {path}: {e}")
        continue

    df["bgc_class"] = os.path.basename(os.path.dirname(path))
    df["source_file"] = path

    frames.append(df)

if not frames:
    raise SystemExit("ERROR: No clustering files loaded.")

all_df = pd.concat(frames, ignore_index=True)

print(f"Total clustering records: {len(all_df)}")

# ------------------------------------------------------------
# Remove MIBiG reference records
# ------------------------------------------------------------

before = len(all_df)

all_df = all_df[
    ~all_df["GBK"].astype(str).str.match(r"^BGC\d{7}$")
].copy()

print(f"MIBiG records removed: {before - len(all_df)}")
print(f"Devosia records retained: {len(all_df)}")

# ------------------------------------------------------------
# Identify GCF column
# ------------------------------------------------------------

possible_gcf = [
    "Family",
    "GCF",
    "family",
    "gcf"
]

gcf_col = None

for c in possible_gcf:
    if c in all_df.columns:
        gcf_col = c
        break

if gcf_col is None:
    raise SystemExit(
        "ERROR: Could not identify GCF column.\n"
        f"Columns: {list(all_df.columns)}"
    )

print(f"Using GCF column: {gcf_col}")

# ------------------------------------------------------------
# Keep selected GCFs
# ------------------------------------------------------------

selected = all_df[
    all_df[gcf_col].astype(str).isin(SELECTED)
].copy()

selected["GCF"] = selected[gcf_col].astype(str)

print("\n=== SELECTED GCF RECORDS ===")
print(
    selected["GCF"]
    .value_counts()
    .sort_index()
    .to_string()
)

# ------------------------------------------------------------
# Parse contig and region
#
# Example:
# NZ_CP186479.1.region001
# ------------------------------------------------------------

def parse_gbk(x):

    x = str(x)

    m = re.match(
        r"^(.+)\.region(\d+)$",
        x
    )

    if m:
        return m.group(1), int(m.group(2))

    return None, None


parsed = selected["GBK"].apply(parse_gbk)

selected["contig"] = parsed.apply(
    lambda x: x[0]
)

selected["region_number"] = parsed.apply(
    lambda x: x[1]
)

print(
    "\nRecords without parsed contig:",
    selected["contig"].isna().sum()
)

# ------------------------------------------------------------
# Load master antiSMASH BGC table
# ------------------------------------------------------------

master_file = "bgc_regions_extended.tsv"

if not os.path.exists(master_file):

    raise SystemExit(
        f"ERROR: Missing {master_file}"
    )

master = pd.read_csv(
    master_file,
    sep="\t"
)

print(
    f"Loaded master BGC table: {len(master)} rows"
)

# ------------------------------------------------------------
# Match using contig + region number
# ------------------------------------------------------------

selected = selected.merge(
    master,
    on=["contig", "region_number"],
    how="left",
    suffixes=("", "_master")
)

# ------------------------------------------------------------
# Reconstruct useful columns
# ------------------------------------------------------------

preferred = [
    "GCF",
    "accession",
    "contig",
    "region_number",
    "start",
    "end",
    "products",
    "candidate_cluster_numbers",
    "contig_edge",
    "bgc_class",
    "GBK",
    "Record",
    "source_file"
]

existing = [
    c for c in preferred
    if c in selected.columns
]

result = selected[existing].copy()

# ------------------------------------------------------------
# Sort
# ------------------------------------------------------------

sort_cols = [
    c for c in
    ["GCF", "accession", "contig", "region_number"]
    if c in result.columns
]

result = result.sort_values(
    sort_cols
)

# ------------------------------------------------------------
# Save combined table
# ------------------------------------------------------------

outfile = (
    "selected_GCF_member_BGCs_c0.7.tsv"
)

result.to_csv(
    outfile,
    sep="\t",
    index=False
)

print(
    f"\nSaved: {outfile}"
)

# ------------------------------------------------------------
# Individual GCF files
# ------------------------------------------------------------

for gcf in sorted(SELECTED):

    sub = result[
        result["GCF"] == gcf
    ].copy()

    out = (
        f"{gcf}_member_BGCs_c0.7.tsv"
    )

    sub.to_csv(
        out,
        sep="\t",
        index=False
    )

    print(
        f"Saved: {out} "
        f"({len(sub)} BGC regions)"
    )

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("\n=== FINAL MEMBER SUMMARY ===")

summary_rows = []

for gcf in sorted(SELECTED):

    sub = result[
        result["GCF"] == gcf
    ]

    n_regions = len(sub)

    n_genomes = (
        sub["accession"]
        .nunique()
        if "accession" in sub.columns
        else None
    )

    classes = (
        ";".join(
            sorted(
                sub["bgc_class"]
                .dropna()
                .astype(str)
                .unique()
            )
        )
        if "bgc_class" in sub.columns
        else ""
    )

    products = (
        ";".join(
            sorted(
                sub["products"]
                .dropna()
                .astype(str)
                .unique()
            )
        )
        if "products" in sub.columns
        else ""
    )

    summary_rows.append({
        "GCF": gcf,
        "BGC_regions": n_regions,
        "Genomes": n_genomes,
        "BGC_classes": classes,
        "Products": products
    })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    "selected_GCF_member_summary_c0.7.tsv",
    sep="\t",
    index=False
)

print(
    summary.to_string(index=False)
)

print(
    "\nSaved: selected_GCF_member_summary_c0.7.tsv"
)

