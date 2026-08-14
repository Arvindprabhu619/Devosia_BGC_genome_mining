
FINAL="$PROJECT_ROOT/final_validation"
TABLE="$FINAL/FINAL_Devosia_Validation.tsv"
FIG="$FINAL/figures"

mkdir -p "$FIG"

echo "============================================================"
echo "GENERATING FINAL DEVOSIA VALIDATION FIGURES"
echo "============================================================"

echo
echo "Input:"
ls -lh "$TABLE"

echo
echo "Output:"
echo "$FIG"

python - "$TABLE" "$FIG" <<'PY'

import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

TABLE = Path(sys.argv[1])
FIG = Path(sys.argv[2])

FIG.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(TABLE, sep="\t")

print("\n============================================================")
print("FINAL DATASET FIGURE GENERATION")
print("============================================================")

print(f"Rows: {len(df)}")
print(f"Columns: {len(df.columns)}")

print("\nColumns:")
for c in df.columns:
    print(" ", c)

# ============================================================
# Convert numeric columns
# ============================================================

numeric_cols = [
    "CheckM2_Completeness",
    "CheckM2_Contamination",
    "Closest_Genome_ANI",
    "Closest_Genome_AF",
    "MSA_Percent",
    "Translation_Table",
    "RED_Value"
]

for c in numeric_cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# ============================================================
# 1. FINAL TAXONOMIC COMPOSITION
# ============================================================

if "GTDB_Genus" in df.columns:

    counts = df["GTDB_Genus"].value_counts()

    plt.figure(figsize=(8, 5))
    counts.sort_values(ascending=True).plot(kind="barh")

    plt.xlabel("Number of genomes")
    plt.ylabel("GTDB-Tk genus")
    plt.title("Final Devosia Dataset: GTDB-Tk Genus Composition")
    plt.tight_layout()

    plt.savefig(
        FIG / "Figure_1_GTDB_genus_distribution.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# ============================================================
# 2. CHECKM2 COMPLETENESS VS CONTAMINATION
# ============================================================

if {
    "CheckM2_Completeness",
    "CheckM2_Contamination"
}.issubset(df.columns):

    plt.figure(figsize=(7, 6))

    plt.scatter(
        df["CheckM2_Completeness"],
        df["CheckM2_Contamination"],
        s=35,
        alpha=0.75
    )

    plt.xlabel("CheckM2 completeness (%)")
    plt.ylabel("CheckM2 contamination (%)")
    plt.title("Final Genomes: CheckM2 Quality")

    plt.axvline(
        90,
        linestyle="--",
        linewidth=1
    )

    plt.axhline(
        5,
        linestyle="--",
        linewidth=1
    )

    plt.tight_layout()

    plt.savefig(
        FIG / "Figure_2_CheckM2_completeness_contamination.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# ============================================================
# 3. COMPLETENESS DISTRIBUTION
# ============================================================

if "CheckM2_Completeness" in df.columns:

    plt.figure(figsize=(7, 5))

    plt.hist(
        df["CheckM2_Completeness"].dropna(),
        bins=15
    )

    plt.xlabel("Completeness (%)")
    plt.ylabel("Number of genomes")
    plt.title("CheckM2 Completeness Distribution")

    plt.tight_layout()

    plt.savefig(
        FIG / "Figure_3_CheckM2_completeness_distribution.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# ============================================================
# 4. CONTAMINATION DISTRIBUTION
# ============================================================

if "CheckM2_Contamination" in df.columns:

    plt.figure(figsize=(7, 5))

    plt.hist(
        df["CheckM2_Contamination"].dropna(),
        bins=15
    )

    plt.xlabel("Contamination (%)")
    plt.ylabel("Number of genomes")
    plt.title("CheckM2 Contamination Distribution")

    plt.tight_layout()

    plt.savefig(
        FIG / "Figure_4_CheckM2_contamination_distribution.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# ============================================================
# 5. SPECIES DISTRIBUTION
# ============================================================

if "GTDB_Species" in df.columns:

    species = df["GTDB_Species"].replace(
        ["", "N/A", "nan"],
        pd.NA
    ).dropna()

    counts = species.value_counts()

    plt.figure(figsize=(10, max(6, len(counts) * 0.22)))

    counts.sort_values(ascending=True).plot(
        kind="barh"
    )

    plt.xlabel("Number of genomes")
    plt.ylabel("GTDB-Tk species")
    plt.title("GTDB-Tk Species Distribution")

    plt.tight_layout()

    plt.savefig(
        FIG / "Figure_5_GTDB_species_distribution.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# ============================================================
# 6. ANI DISTRIBUTION
# ============================================================

if "Closest_Genome_ANI" in df.columns:

    ani = df["Closest_Genome_ANI"].dropna()

    plt.figure(figsize=(7, 5))

    plt.hist(
        ani,
        bins=15
    )

    plt.xlabel("Closest genome ANI (%)")
    plt.ylabel("Number of genomes")
    plt.title("GTDB-Tk Closest Genome ANI Distribution")

    plt.tight_layout()

    plt.savefig(
        FIG / "Figure_6_GTDB_ANI_distribution.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# ============================================================
# 7. ANI VS COMPLETENESS
# ============================================================

if {
    "Closest_Genome_ANI",
    "CheckM2_Completeness"
}.issubset(df.columns):

    plt.figure(figsize=(7, 6))

    plt.scatter(
        df["CheckM2_Completeness"],
        df["Closest_Genome_ANI"],
        s=35,
        alpha=0.75
    )

    plt.xlabel("CheckM2 completeness (%)")
    plt.ylabel("Closest genome ANI (%)")
    plt.title("Genome Quality vs GTDB-Tk ANI")

    plt.tight_layout()

    plt.savefig(
        FIG / "Figure_7_Completeness_vs_ANI.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# ============================================================
# 8. VALIDATION FILTERING SUMMARY
# ============================================================

summary = {
    "Initial genomes": 137,
    "Final retained": 110,
    "CheckM2 excluded": 17,
    "Non-Devosia": 10
}

plt.figure(figsize=(8, 5))

labels = list(summary.keys())
values = list(summary.values())

plt.bar(labels, values)

plt.ylabel("Number of genomes")
plt.title("Devosia Genome Validation Workflow")

plt.xticks(
    rotation=20,
    ha="right"
)

plt.tight_layout()

plt.savefig(
    FIG / "Figure_8_Genome_validation_filtering.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ============================================================
# SAVE SUMMARY TABLE
# ============================================================

summary_df = pd.DataFrame(
    {
        "Metric": [
            "Initial genomes",
            "CheckM2 excluded",
            "Non-Devosia",
            "Final Devosia genomes",
            "GTDB genera",
            "GTDB species assigned"
        ],
        "Value": [
            137,
            17,
            10,
            110,
            df["GTDB_Genus"].nunique()
            if "GTDB_Genus" in df.columns else "NA",
            df["GTDB_Species"].replace(
                ["", "N/A", "nan"],
                pd.NA
            ).dropna().nunique()
            if "GTDB_Species" in df.columns else "NA"
        ]
    }
)

summary_df.to_csv(
    FIG / "Validation_Figure_Summary.tsv",
    sep="\t",
    index=False
)

print("\n============================================================")
print("FIGURES GENERATED")
print("============================================================")

for f in sorted(FIG.glob("*.png")):
    print(f.name)

print("\nSummary:")
print(summary_df.to_string(index=False))

PY

echo
echo "============================================================"
echo "FIGURE FILES"
echo "============================================================"

ls -lh "$FIG"

echo
echo "DONE"
