#!/usr/bin/env python3
# ============================================================
# scripts/07_extract_bgcs.py
# ============================================================
# Purpose: Parse antiSMASH JSON output into structured BGC tables.
#
# Inputs:   results/06_antismash/raw_output/*/  (*.json files)
# Outputs:
#   results/06_antismash/parsed/bgc_regions.tsv
#   results/06_antismash/parsed/bgc_class_assignments.tsv
#   results/06_antismash/parsed/hybrid_regions.tsv
#   results/06_antismash/parsed/summary.txt
#
# Runtime:  5 min
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
    for rec in data.get("records", []):
        contig_id = rec.get("id", "")
        for region in rec.get("areas", []):
            rnum = region.get("region_number", 0)
            start = region.get("start", 0)
            end = region.get("end", 0)
            products = region.get("products", [])
            edge = region.get("contig_edge", False)
            ccs = region.get("candidate_clusters", [])
            n_asdom = sum(len(cc.get("domains", [])) for cc in ccs)
            out.append({
                "bgc_id": f"{acc}_region{rnum:03d}",
                "genome_accession": acc,
                "contig_id": contig_id,
                "region_number": rnum,
                "start": start,
                "end": end,
                "length": end - start,
                "product_classes": ";".join(products),
                "n_products": len(products),
                "is_hybrid": len(products) > 1,
                "contig_edge": edge,
                "n_candidate_clusters": len(ccs),
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
        print("ERROR: No JSON files. Did antiSMASH run?")
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

    # Write outputs
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
    summary = f"""Physical BGC regions: {len(all_regions)}
Hybrid regions:       {len(hybrids)}
Class assignments:    {len(class_assignments)}

Top 10 classes:
"""
    for cls, cnt in class_counts.most_common(10):
        summary += f"  {cls:30s} {cnt}\n"
    print(summary)
    (PARSED / "summary.txt").write_text(summary)


if __name__ == "__main__":
    main()
