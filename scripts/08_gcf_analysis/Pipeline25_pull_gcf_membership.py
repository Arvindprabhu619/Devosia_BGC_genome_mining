#!/usr/bin/env python3
"""
pull_gcf_membership.py (v2)

Builds GCF membership + per-family summary from BiG-SCAPE v2-style output,
using record_annotations.tsv's Organism field as the genome identifier
(GBK values in clustering files are contig+region IDs, not genome accessions).
"""

import sys, re, glob, os, argparse
import pandas as pd

MIBIG_RE = re.compile(r"^BGC\d{7}$")


def find_cutoff_dirs(base_dir):
    parent = os.path.dirname(base_dir.rstrip("/"))
    stem = os.path.basename(base_dir.rstrip("/"))
    candidates = sorted(glob.glob(os.path.join(parent, stem + "_c*")))
    if not candidates:
        candidates = sorted(glob.glob(os.path.join(base_dir, "*_c*")))
    return candidates


def parse_cutoff_from_dirname(d):
    m = re.search(r"_c(\d+\.\d+)$", d)
    return float(m.group(1)) if m else None


def load_clustering_files(cutoff_dir):
    tsvs = glob.glob(os.path.join(cutoff_dir, "*", "*_clustering_c*.tsv"))
    frames = []
    for path in tsvs:
        bgc_class = os.path.basename(os.path.dirname(path))
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception as e:
            print(f"  [WARN] could not read {path}: {e}")
            continue
        df["bgc_class"] = bgc_class
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base_dir")
    ap.add_argument("--n-genomes", type=int, default=110)
    args = ap.parse_args()

    cutoff_dirs = find_cutoff_dirs(args.base_dir)
    if not cutoff_dirs:
        print(f"No cutoff subdirectories found from base: {args.base_dir}")
        sys.exit(1)

    all_membership = []

    for cdir in cutoff_dirs:
        cutoff = parse_cutoff_from_dirname(cdir)
        print(f"=== Processing cutoff {cutoff} ===")

        df = load_clustering_files(cdir)
        if df is None:
            print("  [WARN] no clustering TSVs found, skipping.")
            continue
        print(f"  Loaded {len(df)} rows.")

        # remove MIBiG reference entries
        n_before = len(df)
        df = df[~df["GBK"].astype(str).str.match(MIBIG_RE)].copy()
        n_removed = n_before - len(df)
        if n_removed:
            print(f"  Removed {n_removed} MIBiG reference rows.")

        # load annotation table for this cutoff dir (genome-level metadata)
        ann_path = os.path.join(cdir, "record_annotations.tsv")
        if not os.path.exists(ann_path):
            print(f"  [ERROR] {ann_path} not found -- cannot map genomes. Skipping cutoff.")
            continue
        ann = pd.read_csv(ann_path, sep="\t")

        # merge on GBK to pull in Organism (genome id) + Class/Category/Description
        keep_cols = [c for c in ["GBK", "Organism", "Class", "Category", "Description", "Taxonomy"] if c in ann.columns]
        ann_small = ann[keep_cols].drop_duplicates(subset="GBK")
        df = df.merge(ann_small, on="GBK", how="left")

        n_missing_org = df["Organism"].isna().sum() if "Organism" in df.columns else len(df)
        if n_missing_org:
            print(f"  [WARN] {n_missing_org} rows have no Organism match after merge.")
        n_unique_genomes = df["Organism"].nunique() if "Organism" in df.columns else 0
        print(f"  Unique Organism values in this cutoff's clustering set: {n_unique_genomes} (expect ~{args.n_genomes})")

        df["GCF_ID"] = df["bgc_class"].astype(str) + "_c" + str(cutoff) + "_fam" + df["Family"].astype(str)
        df["cutoff"] = cutoff
        df = df.rename(columns={"Organism": "genome"})

        all_membership.append(df)
        print()

    if not all_membership:
        print("No membership data collected -- aborting.")
        sys.exit(1)

    membership = pd.concat(all_membership, ignore_index=True)
    membership.to_csv("gcf_membership_full.tsv", sep="\t", index=False)
    print(f"Saved: gcf_membership_full.tsv ({len(membership)} rows)")

    # ------------------------------------------------------------------
    # Per-GCF summary
    # ------------------------------------------------------------------
    summary_rows = []
    for (cutoff, gcf_id), grp in membership.groupby(["cutoff", "GCF_ID"]):
        n_genomes_in_fam = grp["genome"].nunique()
        pct = 100 * n_genomes_in_fam / args.n_genomes
        category = "core" if pct >= 90 else ("intermediate" if pct >= 10 else "rare")

        products = None
        if "Description" in grp.columns:
            products = "; ".join(sorted(set(grp["Description"].dropna().astype(str))))

        summary_rows.append({
            "cutoff": cutoff,
            "GCF_ID": gcf_id,
            "bgc_class": grp["bgc_class"].iloc[0],
            "n_BGCs": len(grp),
            "n_genomes": n_genomes_in_fam,
            "pct_genomes": round(pct, 1),
            "category": category,
            "descriptions": products,
            "member_GBKs": "; ".join(grp["GBK"].astype(str).tolist()),
        })

    summary = pd.DataFrame(summary_rows).sort_values(
        ["cutoff", "category", "n_genomes"], ascending=[True, True, False]
    )
    summary.to_csv("gcf_family_summary.tsv", sep="\t", index=False)
    print(f"Saved: gcf_family_summary.tsv ({len(summary)} GCFs)")

    print("\n=== Category counts per cutoff (sanity check vs your earlier table) ===")
    print(summary.groupby(["cutoff", "category"]).size().unstack(fill_value=0))

    print("\n=== Core GCFs found ===")
    core = summary[summary["category"] == "core"]
    if len(core):
        print(core[["cutoff", "GCF_ID", "bgc_class", "n_genomes", "pct_genomes", "descriptions"]].to_string(index=False))
    else:
        print("  (none)")


if __name__ == "__main__":
    main()
