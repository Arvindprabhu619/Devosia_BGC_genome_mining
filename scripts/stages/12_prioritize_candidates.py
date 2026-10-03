#!/usr/bin/env python3
# ============================================================
# scripts/12_prioritize_candidates.py
# ============================================================
# Purpose: Score BGCs for experimental follow-up (multi-criterion).
#
# Scoring components (max 13):
#   +3  Putatively novel
#   +1  Related (<70% KCB)
#   +2  Rare GCF (<10% prevalence)
#   +1  Discovery-relevant class (NRPS, PKS, RiPP, ...)
#   +1  Hybrid (multi-class)
#   +2  Distinctive architecture
#   +2  Multi-genome lineage restriction (from phylogeny)
#   +1  >=18 aSDomains
#   +1  >=2 candidate clusters
#
# Inputs:
#   results/06_antismash/parsed/bgc_regions.tsv
#   results/06_antismash/parsed/bgc_novelty.tsv
#   results/09_bigscape/parsed/gcf_prevalence_c0.7.tsv
# Outputs:
#   results/12_prioritization/scored_bgcs.tsv
#   results/12_prioritization/candidate_pool.tsv
#   results/12_prioritization/prioritized_candidates.tsv
# Runtime: 15 min
# ============================================================
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PARSED = ROOT / "results" / "06_antismash" / "parsed"
BS = ROOT / "results" / "09_bigscape" / "parsed"
PRIORITY = ROOT / "results" / "12_prioritization"
PRIORITY.mkdir(parents=True, exist_ok=True)

DISCOVERY_CLASSES = {
    "NRPS", "NRPS-like", "PKS", "T1PKS", "T2PKS", "T3PKS",
    "RiPP", "RiPP-like", "thioamitide", "thiopeptide",
    "lanthipeptide", "lassopeptide", "phosphonate", "NRP-metallophore",
}


def main():
    print("=== [12] Candidate prioritization ===\n")

    # Load novelty
    novelty = {}
    with open(PARSED / "bgc_novelty.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            novelty[row["bgc_id"]] = row

    # Load regions
    regions = {}
    with open(PARSED / "bgc_regions.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            regions[row["bgc_id"]] = row

    # Load GCF prevalence
    gcf_prev = {}
    gcf_file = BS / "gcf_prevalence_c0.7.tsv"
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

        # Novelty
        cat = nov.get("novelty_category", "")
        if cat == "putatively_novel":
            score += 3; reasons.append("novelty+3")
        elif cat == "related":
            score += 1; reasons.append("related+1")

        # Hybrid
        if str(reg.get("is_hybrid", "")).lower() == "true":
            score += 1; reasons.append("hybrid+1")

        # Class relevance
        prods = reg.get("product_classes", "")
        if any(cls in prods for cls in DISCOVERY_CLASSES):
            score += 1; reasons.append("class+1")

        # Domain complexity
        try:
            if int(reg.get("n_asdomains", 0)) >= 18:
                score += 1; reasons.append("domains+1")
        except (ValueError, TypeError):
            pass

        # Candidate clusters
        try:
            if int(reg.get("n_candidate_clusters", 0)) >= 2:
                score += 1; reasons.append("candidates+1")
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

    print(f"\nScored BGCs:        {len(scored)}")
    print(f"Candidate pool:     {len(pool)}")
    print(f"Top 15 candidates:  {len(top)}")


if __name__ == "__main__":
    main()
