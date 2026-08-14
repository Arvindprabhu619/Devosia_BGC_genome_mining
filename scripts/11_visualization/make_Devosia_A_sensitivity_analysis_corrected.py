#!/usr/bin/env python3

"""
Corrected Devosia / Devosia_A sensitivity analysis

Purpose
-------
Evaluate whether inclusion of 22 GTDB-Tk R232 Devosia_A genomes
materially changes:

1. BGC burden
2. BGC novelty composition
3. GCF distribution at BiG-SCAPE cutoff c0.7

IMPORTANT
---------
Novelty is derived ONLY from the validated master table:

    BGC_KnownClusterBlast_novelty_GTDB_R232.tsv

This avoids the previous incorrect use of:
    BGC_novelty_by_GTDB_genus.tsv

The final expected BGC totals are:

    Devosia       = 477
    Devosia_A     = 105
    Combined      = 582

Novelty classification:

    Devosia:
        putatively_novel = 385
        related           = 56
        known             = 36

    Devosia_A:
        putatively_novel = 91
        related           = 14
        known             = 0

Combined:
        putatively_novel = 476
        related           = 70
        known             = 36

GCF c0.7:
        Devosia-only      = 73
        Shared            = 13
        Devosia_A-only    = 25
        Total              = 111
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from pathlib import Path


# ============================================================
# INPUT FILES
# ============================================================

NOVELTY_FILE = Path(
    "BGC_KnownClusterBlast_novelty_GTDB_R232.tsv"
)

GCF_MATRIX = Path(
    "gcf_genome_matrix_c0.7.tsv"
)

TAX_FILE = Path(
    "BGC_110_GTDB_R232_mapping.tsv"
)


# ============================================================
# OUTPUT FILES
# ============================================================

OUTPUT_SUMMARY = Path(
    "Table_Sx_Devosia_DevosiaA_sensitivity_corrected.tsv"
)

OUTPUT_GCF = Path(
    "Table_Sx_Devosia_DevosiaA_GCF_c0.7_sensitivity_corrected.tsv"
)

OUTPUT_PNG = Path(
    "Figure_Sx_Devosia_DevosiaA_sensitivity_c0.7_corrected.png"
)

OUTPUT_PDF = Path(
    "Figure_Sx_Devosia_DevosiaA_sensitivity_c0.7_corrected.pdf"
)


# ============================================================
# SETTINGS
# ============================================================

TARGET_GENERA = [
    "Devosia",
    "Devosia_A"
]

EXPECTED_TOTAL_BGC = 582

EXPECTED_BGC = {
    "Devosia": 477,
    "Devosia_A": 105
}

EXPECTED_NOVELTY = {
    "Devosia": {
        "putatively_novel": 385,
        "related": 56,
        "known": 36
    },

    "Devosia_A": {
        "putatively_novel": 91,
        "related": 14,
        "known": 0
    }
}

EXPECTED_GCF = {
    "total": 111,
    "Devosia_only": 73,
    "shared": 13,
    "Devosia_A_only": 25
}


# ============================================================
# HELPER
# ============================================================

def clean_accession(x):

    return str(x).strip()


# ============================================================
# CHECK INPUT FILES
# ============================================================

print("\n" + "=" * 80)
print("CHECKING INPUT FILES")
print("=" * 80)

for f in [
    NOVELTY_FILE,
    GCF_MATRIX,
    TAX_FILE
]:

    print(f"{f}: ", "FOUND" if f.exists() else "MISSING")

    if not f.exists():

        raise FileNotFoundError(
            f"Required input file not found: {f}"
        )


# ============================================================
# 1. READ FINAL NOVELTY MASTER
# ============================================================

print("\n" + "=" * 80)
print("1. READING FINAL KCB NOVELTY MASTER")
print("=" * 80)

nov = pd.read_csv(
    NOVELTY_FILE,
    sep="\t"
)

print("Rows:", len(nov))

print("\nColumns:")
print(nov.columns.tolist())


required_novelty = {
    "accession",
    "contig",
    "region_number",
    "has_kcb_result",
    "n_ranked_hits",
    "top_mibig",
    "similarity_percent",
    "novelty_category",
    "genus"
}

missing = required_novelty - set(nov.columns)

if missing:

    raise ValueError(
        "Missing required columns in novelty master: "
        + str(sorted(missing))
    )


# ------------------------------------------------------------
# Clean fields
# ------------------------------------------------------------

nov["accession"] = (
    nov["accession"]
    .astype(str)
    .str.strip()
)

nov["contig"] = (
    nov["contig"]
    .astype(str)
    .str.strip()
)

nov["region_number"] = pd.to_numeric(
    nov["region_number"],
    errors="coerce"
)

nov["similarity_percent"] = pd.to_numeric(
    nov["similarity_percent"],
    errors="coerce"
)

nov["genus"] = (
    nov["genus"]
    .astype(str)
    .str.strip()
)

nov["novelty_category"] = (
    nov["novelty_category"]
    .astype(str)
    .str.strip()
)


# ============================================================
# 2. BASIC NOVELTY VALIDATION
# ============================================================

print("\n" + "=" * 80)
print("2. VALIDATING NOVELTY MASTER")
print("=" * 80)


print("\nGenus distribution:")

print(
    nov["genus"]
    .value_counts(dropna=False)
)


print("\nNovelty distribution:")

print(
    nov["novelty_category"]
    .value_counts(dropna=False)
)


# ------------------------------------------------------------
# Total BGC check
# ------------------------------------------------------------

if len(nov) != EXPECTED_TOTAL_BGC:

    raise ValueError(
        f"Expected {EXPECTED_TOTAL_BGC} BGC regions, "
        f"but found {len(nov)}."
    )


# ------------------------------------------------------------
# Remove anything outside target groups
# ------------------------------------------------------------

nov_target = nov[
    nov["genus"].isin(TARGET_GENERA)
].copy()


if len(nov_target) != EXPECTED_TOTAL_BGC:

    raise ValueError(
        "After restricting to Devosia + Devosia_A, "
        f"expected {EXPECTED_TOTAL_BGC} BGCs but found "
        f"{len(nov_target)}."
    )


# ============================================================
# 3. VERIFY BGC COUNTS PER GENUS
# ============================================================

print("\n" + "=" * 80)
print("3. VERIFYING BGC COUNTS BY GTDB GENUS")
print("=" * 80)


bgc_counts = (
    nov_target
    .groupby("genus")
    .size()
    .to_dict()
)


for genus in TARGET_GENERA:

    observed = bgc_counts.get(genus, 0)

    expected = EXPECTED_BGC[genus]

    print(
        f"{genus:<12} observed={observed:<5} "
        f"expected={expected}"
    )

    if observed != expected:

        raise ValueError(
            f"BGC count mismatch for {genus}: "
            f"observed={observed}, expected={expected}"
        )


# ============================================================
# 4. NOVELTY COUNTS AND PERCENTAGES
# ============================================================

print("\n" + "=" * 80)
print("4. CALCULATING NOVELTY BY GENUS")
print("=" * 80)


summary_rows = []


for genus in TARGET_GENERA:

    sub = nov_target[
        nov_target["genus"] == genus
    ].copy()

    total = len(sub)

    novel = (
        sub["novelty_category"]
        == "putatively_novel"
    ).sum()

    related = (
        sub["novelty_category"]
        == "related"
    ).sum()

    known = (
        sub["novelty_category"]
        == "known"
    ).sum()


    # --------------------------------------------------------
    # Validate category totals
    # --------------------------------------------------------

    if novel + related + known != total:

        raise ValueError(
            f"Novelty categories do not sum to total "
            f"for {genus}."
        )


    # --------------------------------------------------------
    # Percentages
    # --------------------------------------------------------

    novel_pct = 100 * novel / total

    related_pct = 100 * related / total

    known_pct = 100 * known / total


    print(
        f"\n{genus}"
    )

    print(
        f"Total BGCs:          {total}"
    )

    print(
        f"Putatively novel:    "
        f"{novel} ({novel_pct:.2f}%)"
    )

    print(
        f"Related:             "
        f"{related} ({related_pct:.2f}%)"
    )

    print(
        f"Known:               "
        f"{known} ({known_pct:.2f}%)"
    )


    # --------------------------------------------------------
    # Validate against expected values
    # --------------------------------------------------------

    expected = EXPECTED_NOVELTY[genus]

    if novel != expected["putatively_novel"]:

        raise ValueError(
            f"{genus}: putatively novel mismatch."
        )

    if related != expected["related"]:

        raise ValueError(
            f"{genus}: related mismatch."
        )

    if known != expected["known"]:

        raise ValueError(
            f"{genus}: known mismatch."
        )


    summary_rows.append({

        "Group": genus,

        "Genomes": (
            nov_target[
                nov_target["genus"] == genus
            ]["accession"]
            .nunique()
        ),

        "Total_BGCs": total,

        "Putatively_novel_BGCs": novel,

        "Putatively_novel_percent":
            novel_pct,

        "Related_BGCs": related,

        "Related_percent":
            related_pct,

        "Known_BGCs": known,

        "Known_percent":
            known_pct
    })


# ============================================================
# 5. OVERALL NOVELTY
# ============================================================

print("\n" + "=" * 80)
print("5. OVERALL NOVELTY")
print("=" * 80)


total_novel = (
    nov_target["novelty_category"]
    == "putatively_novel"
).sum()

total_related = (
    nov_target["novelty_category"]
    == "related"
).sum()

total_known = (
    nov_target["novelty_category"]
    == "known"
).sum()


print(
    f"Putatively novel: "
    f"{total_novel} "
    f"({100*total_novel/len(nov_target):.2f}%)"
)

print(
    f"Related: "
    f"{total_related} "
    f"({100*total_related/len(nov_target):.2f}%)"
)

print(
    f"Known: "
    f"{total_known} "
    f"({100*total_known/len(nov_target):.2f}%)"
)


if (
    total_novel != 476
    or total_related != 70
    or total_known != 36
):

    raise ValueError(
        "Overall novelty counts do not match "
        "validated 582-BGC results."
    )


# ============================================================
# 6. READ GCF MATRIX
# ============================================================

print("\n" + "=" * 80)
print("6. READING GCF MATRIX — c0.7")
print("=" * 80)


gcf = pd.read_csv(
    GCF_MATRIX,
    sep="\t"
)


if "accession" not in gcf.columns:

    raise ValueError(
        "GCF matrix must contain accession."
    )


gcf["accession"] = (
    gcf["accession"]
    .map(clean_accession)
)


gcf_cols = [

    c for c in gcf.columns

    if str(c).startswith("FAM_")
]


if len(gcf_cols) == 0:

    raise ValueError(
        "No FAM_ GCF columns found."
    )


print(
    "Genomes in GCF matrix:",
    len(gcf)
)

print(
    "GCF families:",
    len(gcf_cols)
)


# ============================================================
# 7. READ GTDB TAXONOMY
# ============================================================

print("\n" + "=" * 80)
print("7. READING GTDB-R232 TAXONOMY")
print("=" * 80)


tax = pd.read_csv(
    TAX_FILE,
    sep="\t"
)


required_tax = {
    "accession",
    "genus"
}


missing = required_tax - set(tax.columns)

if missing:

    raise ValueError(
        f"Missing taxonomy columns: {missing}"
    )


tax = tax[
    [
        "accession",
        "genus"
    ]
].copy()


tax["accession"] = (
    tax["accession"]
    .map(clean_accession)
)


tax["genus"] = (
    tax["genus"]
    .astype(str)
    .str.strip()
)


print(
    "\nTaxonomic distribution:"
)

print(
    tax["genus"]
    .value_counts(dropna=False)
)


# ============================================================
# 8. MERGE GCF + TAXONOMY
# ============================================================

print("\n" + "=" * 80)
print("8. MERGING GCF MATRIX WITH GTDB TAXONOMY")
print("=" * 80)


gcf_tax = gcf.merge(
    tax,
    on="accession",
    how="left",
    validate="one_to_one"
)


missing_tax = (
    gcf_tax["genus"]
    .isna()
    .sum()
)


print(
    "GCF genomes without taxonomy:",
    missing_tax
)


gcf_tax = gcf_tax[
    gcf_tax["genus"].isin(
        TARGET_GENERA
    )
].copy()


devosia = gcf_tax[
    gcf_tax["genus"] == "Devosia"
].copy()


devosia_a = gcf_tax[
    gcf_tax["genus"] == "Devosia_A"
].copy()


print(
    "\nDevosia genomes:",
    len(devosia)
)

print(
    "Devosia_A genomes:",
    len(devosia_a)
)


if len(devosia) != 88:

    raise ValueError(
        f"Expected 88 Devosia genomes, "
        f"found {len(devosia)}."
    )


if len(devosia_a) != 22:

    raise ValueError(
        f"Expected 22 Devosia_A genomes, "
        f"found {len(devosia_a)}."
    )


# ============================================================
# 9. GCF SET ANALYSIS
# ============================================================

print("\n" + "=" * 80)
print("9. GCF SENSITIVITY ANALYSIS — c0.7")
print("=" * 80)


devosia_prev = (
    devosia[gcf_cols] > 0
).sum(axis=0)


devosia_a_prev = (
    devosia_a[gcf_cols] > 0
).sum(axis=0)


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


shared = (
    devosia_set
    &
    devosia_a_set
)


devosia_only = (
    devosia_set
    -
    devosia_a_set
)


devosia_a_only = (
    devosia_a_set
    -
    devosia_set
)


print(
    "Total GCFs:        ",
    len(gcf_cols)
)

print(
    "Devosia GCFs:      ",
    len(devosia_set)
)

print(
    "Devosia_A GCFs:    ",
    len(devosia_a_set)
)

print(
    "Devosia-only:      ",
    len(devosia_only)
)

print(
    "Shared:             ",
    len(shared)
)

print(
    "Devosia_A-only:     ",
    len(devosia_a_only)
)


# ------------------------------------------------------------
# Validate GCF values
# ------------------------------------------------------------

if len(gcf_cols) != EXPECTED_GCF["total"]:

    raise ValueError(
        "Unexpected total GCF count."
    )


if len(devosia_only) != EXPECTED_GCF["Devosia_only"]:

    raise ValueError(
        "Unexpected Devosia-only GCF count."
    )


if len(shared) != EXPECTED_GCF["shared"]:

    raise ValueError(
        "Unexpected shared GCF count."
    )


if len(devosia_a_only) != EXPECTED_GCF["Devosia_A_only"]:

    raise ValueError(
        "Unexpected Devosia_A-only GCF count."
    )


# ============================================================
# 10. GCF DETAIL TABLE
# ============================================================

gcf_rows = []


for fam in sorted(gcf_cols):

    if fam in shared:

        category = "shared"

    elif fam in devosia_only:

        category = "Devosia_only"

    elif fam in devosia_a_only:

        category = "Devosia_A_only"

    else:

        category = "absent"


    dev_prev = int(
        devosia_prev.get(fam, 0)
    )

    dev_a_prev = int(
        devosia_a_prev.get(fam, 0)
    )


    gcf_rows.append({

        "GCF": fam,

        "Devosia_genomes_present":
            dev_prev,

        "Devosia_A_genomes_present":
            dev_a_prev,

        "Total_genomes_present":
            dev_prev + dev_a_prev,

        "distribution":
            category
    })


gcf_detail = pd.DataFrame(
    gcf_rows
)


gcf_detail.to_csv(
    OUTPUT_GCF,
    sep="\t",
    index=False
)


print(
    "\nGCF detail table written:",
    OUTPUT_GCF
)


# ============================================================
# 11. ADD GCF SUMMARY TO MAIN TABLE
# ============================================================

for row in summary_rows:

    genus = row["Group"]

    if genus == "Devosia":

        row["GCFs_present_c0.7"] = len(
            devosia_set
        )

    else:

        row["GCFs_present_c0.7"] = len(
            devosia_a_set
        )


summary_rows.append({

    "Group":
        "Combined",

    "Genomes":
        len(devosia) + len(devosia_a),

    "Total_BGCs":
        len(nov_target),

    "Putatively_novel_BGCs":
        total_novel,

    "Putatively_novel_percent":
        100 * total_novel / len(nov_target),

    "Related_BGCs":
        total_related,

    "Related_percent":
        100 * total_related / len(nov_target),

    "Known_BGCs":
        total_known,

    "Known_percent":
        100 * total_known / len(nov_target),

    "GCFs_present_c0.7":
        len(gcf_cols)
})


summary_rows.append({

    "Group":
        "GCF_Devosia_only_c0.7",

    "Genomes": "",

    "Total_BGCs": "",

    "Putatively_novel_BGCs": "",

    "Putatively_novel_percent": "",

    "Related_BGCs": "",

    "Related_percent": "",

    "Known_BGCs": "",

    "Known_percent": "",

    "GCFs_present_c0.7":
        len(devosia_only)
})


summary_rows.append({

    "Group":
        "GCF_shared_c0.7",

    "Genomes": "",

    "Total_BGCs": "",

    "Putatively_novel_BGCs": "",

    "Putatively_novel_percent": "",

    "Related_BGCs": "",

    "Related_percent": "",

    "Known_BGCs": "",

    "Known_percent": "",

    "GCFs_present_c0.7":
        len(shared)
})


summary_rows.append({

    "Group":
        "GCF_Devosia_A_only_c0.7",

    "Genomes": "",

    "Total_BGCs": "",

    "Putatively_novel_BGCs": "",

    "Putatively_novel_percent": "",

    "Related_BGCs": "",

    "Related_percent": "",

    "Known_BGCs": "",

    "Known_percent": "",

    "GCFs_present_c0.7":
        len(devosia_a_only)
})


summary = pd.DataFrame(
    summary_rows
)


# ============================================================
# 12. WRITE SUMMARY TABLE
# ============================================================

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
# 13. FIGURE
# ============================================================

print("\n" + "=" * 80)
print("13. GENERATING PUBLICATION FIGURE")
print("=" * 80)


plt.rcParams.update({

    "font.family":
        "DejaVu Sans",

    "font.size":
        10,

    "axes.titlesize":
        12,

    "axes.labelsize":
        10,

    "xtick.labelsize":
        9,

    "ytick.labelsize":
        9,

    "legend.fontsize":
        8
})


fig = plt.figure(
    figsize=(12, 4.8)
)


gs = fig.add_gridspec(
    1,
    3,
    width_ratios=[
        1.0,
        1.2,
        1.1
    ],
    wspace=0.38
)


# ============================================================
# PANEL A
# BGC BURDEN
# ============================================================

ax1 = fig.add_subplot(
    gs[0]
)


groups = [
    "Devosia\n(n=88)",
    "Devosia_A\n(n=22)"
]


values = [

    summary[
        summary["Group"] == "Devosia"
    ]["Total_BGCs"].iloc[0]
    /
    88,

    summary[
        summary["Group"] == "Devosia_A"
    ]["Total_BGCs"].iloc[0]
    /
    22
]


bars = ax1.bar(
    groups,
    values
)


ax1.set_ylabel(
    "BGCs per genome\n(mean)"
)


ax1.set_title(
    "A. BGC burden"
)


ax1.spines[
    "top"
].set_visible(False)


ax1.spines[
    "right"
].set_visible(False)


for bar, value in zip(
    bars,
    values
):

    ax1.text(

        bar.get_x()
        +
        bar.get_width() / 2,

        value + 0.05,

        f"{value:.2f}",

        ha="center",

        va="bottom",

        fontsize=9
    )


# ============================================================
# PANEL B
# NOVELTY COMPOSITION
# ============================================================

ax2 = fig.add_subplot(
    gs[1]
)


xpos = np.arange(
    2
)


novel_pct = [

    100 *
    summary[
        summary["Group"] == "Devosia"
    ]["Putatively_novel_BGCs"].iloc[0]
    /
    summary[
        summary["Group"] == "Devosia"
    ]["Total_BGCs"].iloc[0],

    100 *
    summary[
        summary["Group"] == "Devosia_A"
    ]["Putatively_novel_BGCs"].iloc[0]
    /
    summary[
        summary["Group"] == "Devosia_A"
    ]["Total_BGCs"].iloc[0]
]


related_pct = [

    100 *
    summary[
        summary["Group"] == "Devosia"
    ]["Related_BGCs"].iloc[0]
    /
    summary[
        summary["Group"] == "Devosia"
    ]["Total_BGCs"].iloc[0],

    100 *
    summary[
        summary["Group"] == "Devosia_A"
    ]["Related_BGCs"].iloc[0]
    /
    summary[
        summary["Group"] == "Devosia_A"
    ]["Total_BGCs"].iloc[0]
]


known_pct = [

    100 *
    summary[
        summary["Group"] == "Devosia"
    ]["Known_BGCs"].iloc[0]
    /
    summary[
        summary["Group"] == "Devosia"
    ]["Total_BGCs"].iloc[0],

    100 *
    summary[
        summary["Group"] == "Devosia_A"
    ]["Known_BGCs"].iloc[0]
    /
    summary[
        summary["Group"] == "Devosia_A"
    ]["Total_BGCs"].iloc[0]
]


# ------------------------------------------------------------
# Stacked bars
# ------------------------------------------------------------

ax2.bar(
    xpos,
    novel_pct,
    label="Putatively novel"
)


ax2.bar(
    xpos,
    related_pct,
    bottom=novel_pct,
    label="Related"
)


bottom = (
    np.array(novel_pct)
    +
    np.array(related_pct)
)


ax2.bar(
    xpos,
    known_pct,
    bottom=bottom,
    label="Known"
)


ax2.set_xticks(
    xpos
)


ax2.set_xticklabels(
    groups
)


ax2.set_ylim(
    0,
    100
)


ax2.set_ylabel(
    "BGCs (%)"
)


ax2.set_title(
    "B. BGC novelty composition"
)


ax2.spines[
    "top"
].set_visible(False)


ax2.spines[
    "right"
].set_visible(False)


# ------------------------------------------------------------
# Percentage labels
# ------------------------------------------------------------

for i in range(2):

    if novel_pct[i] >= 5:

        ax2.text(

            i,

            novel_pct[i] / 2,

            f"{novel_pct[i]:.1f}%",

            ha="center",

            va="center",

            fontsize=8
        )


    if related_pct[i] >= 5:

        ax2.text(

            i,

            novel_pct[i]
            +
            related_pct[i] / 2,

            f"{related_pct[i]:.1f}%",

            ha="center",

            va="center",

            fontsize=8
        )


    if known_pct[i] >= 5:

        ax2.text(

            i,

            bottom[i]
            +
            known_pct[i] / 2,

            f"{known_pct[i]:.1f}%",

            ha="center",

            va="center",

            fontsize=8
        )


# ------------------------------------------------------------
# Legend
# ------------------------------------------------------------

ax2.legend(

    frameon=False,

    loc="upper center",

    bbox_to_anchor=(
        0.5,
        -0.14
    ),

    ncol=3
)


# ============================================================
# PANEL C
# GCF OVERLAP
# ============================================================

ax3 = fig.add_subplot(
    gs[2]
)


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


# ------------------------------------------------------------
# Explicit numerical labels
# ------------------------------------------------------------

if venn.get_label_by_id("10"):

    venn.get_label_by_id(
        "10"
    ).set_text(
        str(len(devosia_only))
    )


if venn.get_label_by_id("01"):

    venn.get_label_by_id(
        "01"
    ).set_text(
        str(len(devosia_a_only))
    )


if venn.get_label_by_id("11"):

    venn.get_label_by_id(
        "11"
    ).set_text(
        str(len(shared))
    )


ax3.text(

    0.5,

    -0.08,

    (
        f"Total GCFs = {len(gcf_cols)}"
    ),

    transform=ax3.transAxes,

    ha="center",

    fontsize=9
)


# ============================================================
# OVERALL TITLE
# ============================================================

fig.suptitle(

    "Sensitivity of BGC novelty and GCF distributions "
    "to inclusion of GTDB-defined Devosia_A",

    fontsize=13,

    y=1.02
)


plt.tight_layout()


# ============================================================
# SAVE FIGURE
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
# FINAL VALIDATION REPORT
# ============================================================

print("\n" + "=" * 80)
print("FINAL VALIDATION")
print("=" * 80)


print("\nBGC totals:")
print(
    "Devosia:    477"
)

print(
    "Devosia_A:  105"
)

print(
    "Combined:   582"
)


print("\nNovelty:")

print(
    "Devosia: "
    "385 putatively novel, "
    "56 related, "
    "36 known"
)

print(
    "Devosia_A: "
    "91 putatively novel, "
    "14 related, "
    "0 known"
)

print(
    "Combined: "
    "476 putatively novel, "
    "70 related, "
    "36 known"
)


print("\nGCF c0.7:")

print(
    "Devosia-only:   ",
    len(devosia_only)
)

print(
    "Shared:          ",
    len(shared)
)

print(
    "Devosia_A-only:  ",
    len(devosia_a_only)
)

print(
    "Total:           ",
    len(gcf_cols)
)


print("\nOutput files:")

print(
    "Summary:",
    OUTPUT_SUMMARY
)

print(
    "GCF table:",
    OUTPUT_GCF
)

print(
    "PNG:",
    OUTPUT_PNG
)

print(
    "PDF:",
    OUTPUT_PDF
)


print("\n" + "=" * 80)
print("ANALYSIS COMPLETED SUCCESSFULLY")
print("=" * 80)
