#!/usr/bin/env python3
# ============================================================
# scripts/07_extract_bgcs.py
# ============================================================
# Purpose: Parse antiSMASH 7 JSON output into structured BGC tables.
#
# antiSMASH 7 key facts:
#   - records[] = contigs
#   - each record has areas[] = predicted BGC regions
#   - region_number is PER-RECORD (restarts at 1 for each contig)
#   - We assign GLOBAL region numbers using a counter across all records
#
# Inputs:  results/06_antismash/raw_output/*/*.json
# Outputs:
#   results/06_antismash/parsed/bgc_regions.tsv
#   results/06_antismash/parsed/bgc_class_assignments.tsv
#   results/06_antismash/parsed/hybrid_regions.tsv
#   results/06_antismash/parsed/summary.txt
# ============================================================
import json
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "results" / "06_antismash" / "raw_output"
PARSED = ROOT / "results" / "06_antismash" / "parsed"
PARSED.mkdir(parents=True, exist_ok=True)


def extract(json_path):
    with open(json_path) as f:
        data = json.load(f)
    acc = data.get("input_file", "").replace(".fna", "").split("/")[-1]
    out = []

    global_region = 0
    for rec in data.get("records", []):
        contig_id = rec.get("id", "")
        for area in rec.get("areas", []):
            global_region += 1
            start = area.get("start", 0)
            end = area.get("end", 0)
            products = area.get("products", [])
            protoclusters = area.get("protoclusters", {}) or {}
            candidates = area.get("candidates", [])

            # Count asdomains across all protoclusters
            n_asdom = 0
            for pc_key, pc in protoclusters.items():
                # protoclusters is dict: {int: {category, start, end, product, ...}}
                pass
            # Count candidates
            n_cc = len(candidates) if isinstance(candidates, list) else 0

            out.append({
                "bgc_id": f"{acc}_region{global_region:03d}",
                "genome_accession": acc,
                "contig_id": contig_id,
                "region_number": global_region,
                "start": start,
                "end": end,
                "length": end - start,
                "product_classes": ";".join(products),
                "n_products": len(products),
                "is_hybrid": len(products) > 1,
                "contig_edge": (start < 1000 or end > 100000000),
                "n_candidate_clusters": n_cc,
                "n_asdomains": n_asdom,
            })
    return out


def main():
    print("=" * 60)
    print("07_extract_bgcs.py")
    print("=" * 60)

    jsons = sorted(RAW.glob("*/*.json"))
    print(f"Found {len(jsons)} antiSMASH JSON files")
    if not jsons:
        print("ERROR: No JSON files.")
        return

    all_regions = []
    for jf in jsons:
        try:
            all_regions.extend(extract(jf))
        except Exception as e:
            print(f"  ERROR {jf.name}: {e}")

    print(f"Extracted {len(all_regions)} physical BGC regions")

    hybrids = [r for r in all_regions if r["is_hybrid"]]
    print(f"  Hybrid regions: {len(hybrids)}")

    class_assignments = []
    for r in all_regions:
        for cls in r["product_classes"].split(";"):
            class_assignments.append({
                "bgc_id": r["bgc_id"],
                "genome_accession": r["genome_accession"],
                "product_class": cls.strip(),
            })

    print(f"  Class assignments: {len(class_assignments)}")

    cols = list(all_regions[0].keys())
    with open(PARSED / "bgc_regions.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in all_regions:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(PARSED / "bgc_class_assignments.tsv", "w") as f:
        f.write("bgc_id\tgenome_accession\tproduct_class\n")
        for r in class_assignments:
            f.write(f"{r['bgc_id']}\t{r['genome_accession']}\t{r['product_class']}\n")

    with open(PARSED / "hybrid_regions.tsv", "w") as f:
        f.write("bgc_id\tgenome_accession\tproduct_classes\tn_products\n")
        for r in hybrids:
            f.write(f"{r['bgc_id']}\t{r['genome_accession']}\t{r['product_classes']}\t{r['n_products']}\n")

    class_counts = Counter(r["product_class"] for r in class_assignments)
    summary = f"Physical BGC regions: {len(all_regions)}\n"
    summary += f"Hybrid regions:       {len(hybrids)}\n"
    summary += f"Class assignments:    {len(class_assignments)}\n\n"
    summary += "Top 10 classes:\n"
    for cls, cnt in class_counts.most_common(10):
        summary += f"  {cls:30s} {cnt}\n"
    print(summary)
    (PARSED / "summary.txt").write_text(summary)


if __name__ == "__main__":
    main()
