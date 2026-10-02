#!/usr/bin/env python3
"""12_prioritize_candidates.py — Score and rank BGCs."""
from pathlib import Path
from collections import defaultdict
import csv

ROOT = Path(__file__).resolve().parent.parent
PARSED = ROOT / "results" / "06_antismash" / "parsed"
BIGSCAPE = ROOT / "results" / "07_bigscape" / "parsed"
PRIORITY = ROOT / "results" / "09_prioritization"
PRIORITY.mkdir(parents=True, exist_ok=True)

S_NOVEL = 3
S_RELATED = 1
S_RARE_GCF = 2
S_DISCOVERY_CLASS = 1
S_HYBRID = 1
S_DOMAIN = 1
S_CANDIDATE = 1

DISCOVERY_CLASSES = {
    "NRPS", "NRPS-like", "PKS", "T1PKS", "T2PKS", "T3PKS",
    "RiPP", "RiPP-like", "thioamitide", "thiopeptide",
    "lanthipeptide", "lassopeptide", "phosphonate", "NRP-metallophore",
}


def main():
    print("=== [12] Candidate prioritization ===\n")

    novelty = {}
    with open(PARSED / "bgc_novelty.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            novelty[row["bgc_id"]] = row

    regions = {}
    with open(PARSED / "bgc_regions.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            regions[row["bgc_id"]] = row

    gcf_prev = {}
    gcf_file = BIGSCAPE / "gcf_prevalence_c0.7.tsv"
    if gcf_file.exists():
        with open(gcf_file) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                gcf_prev[row["gcf_id"]] = row

    print(f"Loaded: {len(novelty)} novelty, {len(regions)} regions, {len(gcf_prev)} GCFs")

    scored = []
    for bgc_id, reg in regions.items():
        nov = novelty.get(bgc_id, {})
        score = 0
        reasons = []

        cat = nov.get("novelty_category", "")
        if cat == "putatively_novel":
            score += S_NOVEL; reasons.append("novelty+3")
        elif cat == "related":
            score += S_RELATED; reasons.append("related+1")

        if str(reg.get("is_hybrid", "")).lower() == "true":
            score += S_HYBRID; reasons.append("hybrid+1")

        prods = reg.get("product_classes", "")
        if any(cls in prods for cls in DISCOVERY_CLASSES):
            score += S_DISCOVERY_CLASS; reasons.append("class+1")

        try:
            if int(reg.get("n_asdomains", 0)) >= 18:
                score += S_DOMAIN; reasons.append("domains+1")
        except (ValueError, TypeError):
            pass

        try:
            if int(reg.get("n_candidate_clusters", 0)) >= 2:
                score += S_CANDIDATE; reasons.append("candidate+1")
        except (ValueError, TypeError):
            pass

        scored.append({
            "bgc_id": bgc_id,
            "genome_accession": reg.get("genome_accession", ""),
            "product_classes": prods,
            "novelty_category": cat,
            "score": score,
            "reasons": ";".join(reasons),
        })

    scored.sort(key=lambda r: -r["score"])

    with open(PRIORITY / "scored_bgcs.tsv", "w") as f:
        f.write("bgc_id\tgenome_accession\tproduct_classes\tnovelty_category\tscore\treasons\n")
        for r in scored:
            f.write(f"{r['bgc_id']}\t{r['genome_accession']}\t{r['product_classes']}\t{r['novelty_category']}\t{r['score']}\t{r['reasons']}\n")

    pool = [r for r in scored if r["score"] >= 6]
    with open(PRIORITY / "candidate_pool.tsv", "w") as f:
        f.write("bgc_id\tgenome_accession\tproduct_classes\tnovelty_category\tscore\treasons\n")
        for r in pool:
            f.write(f"{r['bgc_id']}\t{r['genome_accession']}\t{r['product_classes']}\t{r['novelty_category']}\t{r['score']}\t{r['reasons']}\n")

    top = pool[:15]
    with open(PRIORITY / "prioritized_candidates.tsv", "w") as f:
        f.write("rank\tbgc_id\tgenome_accession\tproduct_classes\tnovelty_category\tscore\treasons\n")
        for i, r in enumerate(top, 1):
            f.write(f"{i}\t{r['bgc_id']}\t{r['genome_accession']}\t{r['product_classes']}\t{r['novelty_category']}\t{r['score']}\t{r['reasons']}\n")

    print(f"Scored BGCs:        {len(scored)}")
    print(f"Candidate pool:     {len(pool)}")
    print(f"Top 15 candidates:  {len(top)}")


if __name__ == "__main__":
    main()
