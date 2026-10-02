#!/usr/bin/env python3
# ============================================================
# scripts/08_classify_novelty.py
# ============================================================
# Purpose: Classify each BGC by KnownClusterBlast similarity to MIBiG.
#
# Inputs:   results/06_antismash/raw_output/*/  (*.json files)
#           results/06_antismash/parsed/bgc_regions.tsv
# Outputs:
#   results/06_antismash/parsed/bgc_novelty.tsv
#   results/06_antismash/parsed/novelty_summary.txt
#
# Classification:
#   Putatively novel: no detectable ranked MIBiG match
#   Related:          match with similarity < 70%
#   Known:            match with similarity >= 70% (study-specific cutoff)
#
# Runtime:  5 min
# ============================================================
import json
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "results" / "06_antismash" / "raw_output"
PARSED = ROOT / "results" / "06_antismash" / "parsed"
SIM_THRESHOLD = 70.0


def get_kcb(json_path):
    with open(json_path) as f:
        data = json.load(f)
    acc = data.get("input_file", "").replace(".fna", "").split("/")[-1]
    out = {}
    for rec in data.get("records", []):
        for region in rec.get("areas", []):
            rnum = region.get("region_number", 0)
            bgc_id = f"{acc}_region{rnum:03d}"
            kcb = region.get("knownclusterblast", {}) or {}
            hits = kcb.get("results", []) or []
            max_sim = 0.0
            best_acc = ""
            if hits:
                max_sim = float(hits[0].get("similarity", 0))
                best_acc = hits[0].get("accession", "")
            out[bgc_id] = {
                "max_similarity": max_sim,
                "mibig_accession": best_acc,
                "n_hits": len(hits),
            }
    return out


def main():
    print("=" * 60)
    print("08_classify_novelty.py")
    print("=" * 60)

    regions = []
    with open(PARSED / "bgc_regions.tsv") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            regions.append(dict(zip(header, line.strip().split("\t"))))
    print(f"Loaded {len(regions)} BGC regions")

    kcb_all = {}
    for jf in RAW.glob("*/*.json"):
        try:
            kcb_all.update(get_kcb(jf))
        except Exception as e:
            print(f"  ERROR {jf.name}: {e}")
    print(f"Found KnownClusterBlast data for {len(kcb_all)} regions")

    output = []
    for r in regions:
        bgc_id = r["bgc_id"]
        kcb = kcb_all.get(bgc_id, {"max_similarity": 0, "mibig_accession": "", "n_hits": 0})
        sim = kcb["max_similarity"]

        if kcb["n_hits"] == 0:
            cat = "putatively_novel"
        elif sim >= SIM_THRESHOLD:
            cat = "known"
        else:
            cat = "related"

        output.append({
            "bgc_id": bgc_id,
            "genome_accession": r["genome_accession"],
            "product_classes": r["product_classes"],
            "is_hybrid": r["is_hybrid"],
            "contig_edge": r["contig_edge"],
            "max_kcb_similarity": sim,
            "mibig_accession": kcb["mibig_accession"],
            "n_kcb_hits": kcb["n_hits"],
            "novelty_category": cat,
        })

    cols = list(output[0].keys())
    with open(PARSED / "bgc_novelty.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in output:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    counts = Counter(r["novelty_category"] for r in output)
    total = len(output)
    summary = f"""Novelty classification (threshold: >= {SIM_THRESHOLD}% = known)

  Putatively novel: {counts.get('putatively_novel', 0):>4d}  ({100*counts.get('putatively_novel', 0)/total:.2f}%)
  Related:          {counts.get('related', 0):>4d}  ({100*counts.get('related', 0)/total:.2f}%)
  Known:            {counts.get('known', 0):>4d}  ({100*counts.get('known', 0)/total:.2f}%)
"""
    print(summary)
    (PARSED / "novelty_summary.txt").write_text(summary)


if __name__ == "__main__":
    main()
