#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


# ============================================================
# INPUT FILES
# ============================================================

NOVELTY_FILE = "BGC_novelty_by_GTDB_genus.tsv"
GCF_MATRIX = "gcf_genome_matrix_c0.7.tsv"
TAX_FILE = "BGC_110_GTDB_R232_mapping.tsv"
BGC_COUNT_FILE = "bgc_counts_per_genome.tsv"

# Output
OUTPUT_PNG = "Figure_Sx_Devosia_DevosiaA_sensitivity_c0.7.png"
OUTPUT_PDF = "Figure_Sx_Devosia_DevosiaA_sensitivity_c0.7.pdf"
OUTPUT_SUMMARY = "Table_Sx_Devosia_DevosiaA_sensitivity_summary.tsv"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def clean_accession(x):
    """
    Clean accession strings while preserving version numbers.
    """
    return str(x).strip()


def read_bgc_count_file(filename):
    """
    Robust reader for BGC count files.

    Handles both:
      1. Headered files:
         accession    BGC_count
         GCF_xxx      5

      2. Headerless files:
         GCF_xxx      5
    """

    # First try normal header
    df = pd.read_csv(filename, sep="\t", dtype=str)

    print("\nRaw BGC count columns:")
    print(list(df.columns))

    # --------------------------------------------------------
    # CASE 1: Header is actually a genome accession
    # --------------------------------------------------------

    first_col = str(df.columns[0])

    if first_col.startswith(("GCF_", "GCA_")):

        print("Detected headerless BGC-count file.")

        # Re-read without header
        df = pd.read_csv(
            filename,
            sep="\t",
            header=None,
            dtype=str
        )

        if df.shape[1] < 2:
            raise ValueError(
                "BGC count file must contain at least two columns: "
                "accession and BGC count."
            )

        df = df.iloc[:, :2].copy()
        df.columns = ["accession", "BGC_count"]

    else:

        print("Detected headered BGC-count file.")

        # ----------------------------------------------------
        # Find accession column
        # ----------------------------------------------------

        accession_col = None

        for c in df.columns:

            values = df[c].astype(str)

            if values.str.startswith(("GCF_", "GCA_")).any():
                accession_col = c
                break

        if accession_col is None:

            # Fallback: assume first column
            accession_col = df.columns[0]

        # ----------------------------------------------------
        # Find numeric BGC count column
        # ----------------------------------------------------

        numeric_candidates = []

        for c in df.columns:

            if c == accession_col:
                continue

            numeric = pd.to_numeric(df[c], errors="coerce")

            if numeric.notna().sum() > 0:
                numeric_candidates.append(c)

        if not numeric_candidates:

            raise ValueError(
                "Could not identify a numeric BGC-count column."
            )

        # Prefer obvious names
        preferred = [
            c for c in numeric_candidates
            if any(
                key in str(c).lower()
                for key in ["bgc", "count", "regions"]
            )
        ]

        if preferred:
            count_col = preferred[0]
        else:
            count_col = numeric_candidates[0]

        df = df[[accession_col, count_col]].copy()

        df.columns = ["accession", "BGC_count"]

    # --------------------------------------------------------
    # Clean
    # --------------------------------------------------------

    df["accession"] = df["accession"].map(clean_accession)

    df["BGC_count"] = pd.to_numeric(
        df["BGC_count"],
        errors="coerce"
    )

    df = df.dropna(subset=["accession", "BGC_count"])

    df["BGC_count"] = df["BGC_count"].astype(float)

    print("\nProcessed BGC count table:")
    print(df.head())

    print("\nNumber of genomes in BGC count table:", len(df))

    return df


# ============================================================
# 1. READ NOVELTY RESULTS
# ============================================================

print("\n" + "=" * 80)
print("READING NOVELTY RESULTS")
print("=" * 80)

nov = pd.read_csv(
    NOVELTY_FILE,
    sep="\t"
)

print(nov.to_string(index=False))

required_novelty = {
    "genus",
    "known",
    "novel",
    "related",
    "Total"
}

missing = required_novelty - set(nov.columns)

if missing:
    raise ValueError(
        f"Missing novelty columns: {missing}"
    )

# Keep only target groups
nov = nov[
    nov["genus"].isin(["Devosia", "Devosia_A"])
].copy()

nov = nov.set_index("genus")


# ============================================================
# 2. READ GCF MATRIX
# ============================================================

print("\n" + "=" * 80)
print("READING GCF MATRIX")
print("=" * 80)

gcf = pd.read_csv(
    GCF_MATRIX,
    sep="\t"
)

if "accession" not in gcf.columns:

    raise ValueError(
        "The GCF matrix must contain an 'accession' column."
    )

gcf["accession"] = gcf["accession"].map(clean_accession)

gcf_cols = [
    c for c in gcf.columns
    if str(c).startswith("FAM_")
]

if len(gcf_cols) == 0:

    raise ValueError(
        "No GCF columns beginning with 'FAM_' were found."
    )

print("Number of genomes:", len(gcf))
print("Number of GCFs:", len(gcf_cols))


# ============================================================
# 3. READ GTDB TAXONOMY
# ============================================================

print("\n" + "=" * 80)
print("READING GTDB-R232 TAXONOMY")
print("=" * 80)

tax = pd.read_csv(
    TAX_FILE,
    sep="\t"
)

required_tax = {"accession", "genus"}

missing = required_tax - set(tax.columns)

if missing:
    raise ValueError(
        f"Missing taxonomy columns: {missing}"
    )

tax = tax[[
    "accession",
    "genus"
]].copy()

tax["accession"] = tax["accession"].map(clean_accession)

print("\nTaxonomic groups:")
print(
    tax["genus"]
    .value_counts(dropna=False)
)


# ============================================================
# 4. MERGE GCF + TAXONOMY
# ============================================================

x = gcf.merge(
    tax,
    on="accession",
    how="left"
)

missing_tax = x["genus"].isna().sum()

if missing_tax > 0:

    print(
        f"\nWARNING: {missing_tax} genomes "
        "have no genus mapping."
    )

x = x[
    x["genus"].isin([
        "Devosia",
        "Devosia_A"
    ])
].copy()

devosia = x[
    x["genus"] == "Devosia"
].copy()

devosia_a = x[
    x["genus"] == "Devosia_A"
].copy()

print("\nFinal sensitivity groups:")
print("Devosia:", len(devosia))
print("Devosia_A:", len(devosia_a))


# ============================================================
# 5. GCF SETS
# ============================================================

devosia_prev = (
    (devosia[gcf_cols] > 0)
    .sum(axis=0)
)

devosia_a_prev = (
    (devosia_a[gcf_cols] > 0)
    .sum(axis=0)
)

devosia_set = set(
    devosia_prev[
        devosia_prev > 0
    ].index
)

devosia_a_set = set(
    devosia_a_prev[
        devosia_a_prev > 0
    ].index
)

shared = devosia_set & devosia_a_set

devosia_only = devosia_set - devosia_a_set

devosia_a_only = devosia_a_set - devosia_set

print("\n" + "=" * 80)
print("GCF SENSITIVITY ANALYSIS — c0.7")
print("=" * 80)

print("Total GCFs:", len(gcf_cols))
print("Devosia GCFs:", len(devosia_set))
print("Devosia_A GCFs:", len(devosia_a_set))
print("Shared:", len(shared))
print("Devosia-only:", len(devosia_only))
print("Devosia_A-only:", len(devosia_a_only))


# ============================================================
# 6. READ BGC COUNTS
# ============================================================

print("\n" + "=" * 80)
print("READING BGC COUNTS")
print("=" * 80)

bgc_counts = read_bgc_count_file(
    BGC_COUNT_FILE
)


# ============================================================
# 7. MERGE BGC COUNTS WITH TAXONOMY
# ============================================================

bgc = bgc_counts.merge(
    tax,
    on="accession",
    how="left"
)

print("\nBGC-count taxonomy distribution:")

print(
    bgc["genus"]
    .value_counts(dropna=False)
)

# Keep only target groups
bgc = bgc[
    bgc["genus"].isin([
        "Devosia",
        "Devosia_A"
    ])
].copy()


# ============================================================
# 8. BGC BURDEN
# ============================================================

bgc_summary = {}

for genus in ["Devosia", "Devosia_A"]:

    sub = bgc[
        bgc["genus"] == genus
    ]

    n = len(sub)

    total_bgcs = sub["BGC_count"].sum()

    mean_bgcs = (
        sub["BGC_count"].mean()
        if n > 0
        else np.nan
    )

    sd_bgcs = (
        sub["BGC_count"].std(ddof=1)
        if n > 1
        else np.nan
    )

    bgc_summary[genus] = {
        "genomes": n,
        "total_bgcs": total_bgcs,
        "mean_bgcs": mean_bgcs,
        "sd_bgcs": sd_bgcs
    }

    print(
        f"{genus}: "
        f"{total_bgcs:.0f} BGCs, "
        f"{mean_bgcs:.2f} ± {sd_bgcs:.2f} BGC/genome"
    )


# ============================================================
# 9. NOVELTY SUMMARY
# ============================================================

novelty_summary = {}

for genus in ["Devosia", "Devosia_A"]:

    row = nov.loc[genus]

    total = float(row["Total"])
    novel = float(row["novel"])
    related = float(row["related"])
    known = float(row["known"])

    novelty_summary[genus] = {
        "total": total,
        "novel": novel,
        "related": related,
        "known": known,
        "novel_pct": novel / total * 100,
        "related_pct": related / total * 100,
        "known_pct": known / total * 100
    }

    print(
        f"\n{genus} novelty:"
    )

    print(
        f"Novel:   {novel:.0f} "
        f"({novel/total*100:.2f}%)"
    )

    print(
        f"Related: {related:.0f} "
        f"({related/total*100:.2f}%)"
    )

    print(
        f"Known:   {known:.0f} "
        f"({known/total*100:.2f}%)"
    )


# ============================================================
# 10. CREATE SUMMARY TABLE
# ============================================================

summary_rows = []

for genus in ["Devosia", "Devosia_A"]:

    summary_rows.append({
        "Group": genus,
        "Genomes": bgc_summary[genus]["genomes"],
        "Total_BGCs": bgc_summary[genus]["total_bgcs"],
        "Mean_BGCs_per_genome":
            bgc_summary[genus]["mean_bgcs"],
        "SD_BGCs_per_genome":
            bgc_summary[genus]["sd_bgcs"],
        "Novel_BGCs":
            novelty_summary[genus]["novel"],
        "Novel_percent":
            novelty_summary[genus]["novel_pct"],
        "Related_BGCs":
            novelty_summary[genus]["related"],
        "Related_percent":
            novelty_summary[genus]["related_pct"],
        "Known_BGCs":
            novelty_summary[genus]["known"],
        "Known_percent":
            novelty_summary[genus]["known_pct"],
        "GCFs_present":
            len(
                devosia_set
                if genus == "Devosia"
                else devosia_a_set
            )
    })

summary_rows.append({
    "Group": "GCF_overlap_c0.7",
    "Genomes": "",
    "Total_BGCs": "",
    "Mean_BGCs_per_genome": "",
    "SD_BGCs_per_genome": "",
    "Novel_BGCs": "",
    "Novel_percent": "",
    "Related_BGCs": "",
    "Related_percent": "",
    "Known_BGCs": "",
    "Known_percent": "",
    "GCFs_present": len(gcf_cols)
})

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUTPUT_SUMMARY,
    sep="\t",
    index=False
)

print(
    "\nSummary table written:",
    OUTPUT_SUMMARY
)


# ============================================================
# 11. FIGURE SETUP
# ============================================================

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "font.family": "DejaVu Sans"
})

fig = plt.figure(
    figsize=(12, 4.5)
)

gs = fig.add_gridspec(
    1,
    3,
    width_ratios=[1, 1.15, 1],
    wspace=0.38
)


# ============================================================
# PANEL A — BGC BURDEN
# ============================================================

ax1 = fig.add_subplot(gs[0])

groups = [
    "Devosia\n(n=88)",
    "Devosia_A\n(n=22)"
]

values = [
    bgc_summary["Devosia"]["mean_bgcs"],
    bgc_summary["Devosia_A"]["mean_bgcs"]
]

errors = [
    bgc_summary["Devosia"]["sd_bgcs"],
    bgc_summary["Devosia_A"]["sd_bgcs"]
]

bars = ax1.bar(
    groups,
    values,
    yerr=errors,
    capsize=4
)

ax1.set_ylabel(
    "BGCs per genome\n(mean ± SD)"
)

ax1.set_title(
    "A. BGC burden"
)

ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

for bar, value in zip(bars, values):

    ax1.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.15,
        f"{value:.2f}",
        ha="center",
        va="bottom",
        fontsize=9
    )


# ============================================================
# PANEL B — NOVELTY
# ============================================================

ax2 = fig.add_subplot(gs[1])

groups = [
    "Devosia\n(n=88)",
    "Devosia_A\n(n=22)"
]

novel_pct = [
    novelty_summary["Devosia"]["novel_pct"],
    novelty_summary["Devosia_A"]["novel_pct"]
]

related_pct = [
    novelty_summary["Devosia"]["related_pct"],
    novelty_summary["Devosia_A"]["related_pct"]
]

known_pct = [
    novelty_summary["Devosia"]["known_pct"],
    novelty_summary["Devosia_A"]["known_pct"]
]

xpos = np.arange(2)

ax2.bar(
    xpos,
    novel_pct,
    label="Novel"
)

ax2.bar(
    xpos,
    related_pct,
    bottom=novel_pct,
    label="Related"
)

bottom = np.array(novel_pct) + np.array(related_pct)

ax2.bar(
    xpos,
    known_pct,
    bottom=bottom,
    label="Known"
)

ax2.set_xticks(xpos)
ax2.set_xticklabels(groups)

ax2.set_ylim(0, 100)

ax2.set_ylabel(
    "BGCs (%)"
)

ax2.set_title(
    "B. BGC novelty"
)

ax2.legend(
    frameon=False,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.12),
    ncol=3
)

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

for i in range(2):

    ax2.text(
        i,
        novel_pct[i] / 2,
        f"{novel_pct[i]:.1f}%",
        ha="center",
        va="center",
        fontsize=8
    )

    ax2.text(
        i,
        novel_pct[i] + related_pct[i] / 2,
        f"{related_pct[i]:.1f}%",
        ha="center",
        va="center",
        fontsize=8
    )

    ax2.text(
        i,
        bottom[i] + known_pct[i] / 2,
        f"{known_pct[i]:.1f}%",
        ha="center",
        va="center",
        fontsize=8
    )


# ============================================================
# PANEL C — GCF OVERLAP
# ============================================================

ax3 = fig.add_subplot(gs[2])

venn = venn2(
    subsets=(
        len(devosia_only),
        len(devosia_a_only),
        len(shared)
    ),
    set_labels=(
        "Devosia",
        "Devosia_A"
    ),
    ax=ax3
)

ax3.set_title(
    "C. GCF overlap at c0.7"
)

# Improve labels
if venn.get_label_by_id("10"):
    venn.get_label_by_id("10").set_text(
        str(len(devosia_only))
    )

if venn.get_label_by_id("01"):
    venn.get_label_by_id("01").set_text(
        str(len(devosia_a_only))
    )

if venn.get_label_by_id("11"):
    venn.get_label_by_id("11").set_text(
        str(len(shared))
    )

ax3.text(
    0.5,
    -0.05,
    f"Total GCFs = {len(gcf_cols)}",
    transform=ax3.transAxes,
    ha="center",
    fontsize=9
)


# ============================================================
# OVERALL FIGURE TITLE
# ============================================================

fig.suptitle(
    "Sensitivity analysis of Devosia_A inclusion",
    fontsize=14,
    y=1.02
)

plt.tight_layout()


# ============================================================
# SAVE
# ============================================================

fig.savefig(
    OUTPUT_PNG,
    dpi=600,
    bbox_inches="tight"
)

fig.savefig(
    OUTPUT_PDF,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# FINAL REPORT
# ============================================================

print("\n" + "=" * 80)
print("FIGURE GENERATED SUCCESSFULLY")
print("=" * 80)

print("PNG:", OUTPUT_PNG)
print("PDF:", OUTPUT_PDF)
print("Summary:", OUTPUT_SUMMARY)

print("\nGCF summary:")
print(
    "Devosia-only :",
    len(devosia_only)
)

print(
    "Shared       :",
    len(shared)
)

print(
    "Devosia_A-only:",
    len(devosia_a_only)
)

print(
    "Total        :",
    len(gcf_cols)
)

print("\nBGC burden:")
for genus in ["Devosia", "Devosia_A"]:

    print(
        f"{genus}: "
        f"{bgc_summary[genus]['total_bgcs']:.0f} BGCs; "
        f"{bgc_summary[genus]['mean_bgcs']:.2f} ± "
        f"{bgc_summary[genus]['sd_bgcs']:.2f} BGC/genome"
    )

print("\nNovelty:")
for genus in ["Devosia", "Devosia_A"]:

    print(
        f"{genus}: "
        f"{novelty_summary[genus]['novel']:.0f} novel "
        f"({novelty_summary[genus]['novel_pct']:.2f}%), "
        f"{novelty_summary[genus]['related']:.0f} related "
        f"({novelty_summary[genus]['related_pct']:.2f}%), "
        f"{novelty_summary[genus]['known']:.0f} known "
        f"({novelty_summary[genus]['known_pct']:.2f}%)"
    )

print("=" * 80)
