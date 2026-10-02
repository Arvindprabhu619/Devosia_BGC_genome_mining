#!/usr/bin/env python3
# ============================================================
# scripts/14_ecological_association.py
# ============================================================
# Purpose: Test BGC class and GCF associations with isolation niches.
#
# Inputs:
#   data/metadata/devosia_refseq.jsonl
#   data/metadata/devosia_genbank.jsonl
#   results/05_gtdbtk/final_genome_list.txt
#   results/06_antismash/parsed/bgc_class_assignments.tsv
#   results/09_bigscape/parsed/gcf_prevalence_c0.7.tsv
# Outputs:
#   results/11_stats/ecological_metadata.tsv
#   results/11_stats/ecological_association_classes.tsv
#   results/11_stats/ecological_association_gcfs.tsv
# Runtime: 30 min
# ============================================================
import json
import re
import csv
from pathlib import Path
from collections import defaultdict

try:
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    HAVE_STATS = True
except ImportError:
    HAVE_STATS = False

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata"
RES_STATS = ROOT / "results" / "11_stats"
RES_STATS.mkdir(parents=True, exist_ok=True)

NICHE_PATTERNS = {
    "soil":            [r"\bsoil\b", r"\brhizosphere\b", r"\bterrestrial\b"],
    "plant":           [r"\bplant\b", r"\brhizoplane\b", r"\bendophyt", r"\bleaf\b", r"\broot\b"],
    "aquatic":         [r"\bmarine\b", r"\bfreshwater\b", r"\bsea\b", r"\bocean\b", r"\briver\b", r"\blake\b"],
    "contaminated":    [r"\bcontaminat", r"\bpollut", r"\bsludge\b", r"\bwaste\b"],
    "host_associated": [r"\bhost\b", r"\bsymbion", r"\bgut\b", r"\binsect\b", r"\bnodule\b"],
    "industrial":      [r"\bindustrial\b", r"\bbioreactor\b", r"\bferment"],
}


def categorize(text):
    if not text:
        return []
    t = text.lower()
    return [n for n, ps in NICHE_PATTERNS.items() if any(re.search(p, t) for p in ps)]


def extract_niche_map():
    m = {}
    for fname in ["devosia_refseq.jsonl", "devosia_genbank.jsonl"]:
        fp = META / fname
        if not fp.exists():
            continue
        with open(fp) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                acc = rec.get("accession", "")
                if not acc:
                    continue
                bs = rec.get("biosample", {}) or {}
                attrs = bs.get("attributes", []) or []
                texts = [a.get("value", "") for a in attrs if a.get("value")]
                ai = rec.get("assembly_info", {}) or {}
                texts.append(ai.get("assembly_name", ""))
                niches = set(categorize(" ".join(texts)))
                if niches:
                    m[acc] = niches
    return m


def test_assoc(genomes, fmap, niches_by, all_niches):
    res = []
    for feat, gset in fmap.items():
        for n in all_niches:
            a = b = c = d = 0
            for g in genomes:
                hf = g in gset
                hn = n in niches_by.get(g, set())
                if hf and hn: a += 1
                elif hf: b += 1
                elif hn: c += 1
                else: d += 1
            if a + c < 3:
                continue
            try:
                odds, p = fisher_exact([[a, b], [c, d]])
                res.append({
                    "feature": feat, "niche": n,
                    "a": a, "b": b, "c": c, "d": d,
                    "odds_ratio": odds, "p_value": p,
                })
            except Exception:
                continue
    return res


def main():
    print("=" * 60)
    print("14_ecological_association.py")
    print("=" * 60)

    final = ROOT / "results" / "05_gtdbtk" / "final_genome_list.txt"
    if not final.exists():
        print(f"ERROR: {final} not found")
        return
    genomes = [g.strip() for g in final.read_text().splitlines() if g.strip()]
    print(f"Genomes: {len(genomes)}")

    niches_map = extract_niche_map()
    niches_by = {g: niches_map.get(g, set()) for g in genomes}
    n_with = sum(1 for v in niches_by.values() if v)
    print(f"With niche info: {n_with}")

    with open(RES_STATS / "ecological_metadata.tsv", "w") as f:
        f.write("accession\tniches\n")
        for g in genomes:
            f.write(f"{g}\t{','.join(sorted(niches_by[g]))}\n")

    if not HAVE_STATS:
        print("scipy/statsmodels not available — metadata saved only")
        return

    # Load class and GCF data
    cls_map = defaultdict(set)
    cf = ROOT / "results" / "06_antismash" / "parsed" / "bgc_class_assignments.tsv"
    if cf.exists():
        with open(cf) as f:
            f.readline()
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 3:
                    cls_map[p[2]].add(p[1])

    gcf_map = defaultdict(set)
    gf = ROOT / "results" / "09_bigscape" / "parsed" / "gcf_prevalence_c0.7.tsv"
    if gf.exists():
        with open(gf) as f:
            h = f.readline().strip().split("\t")
            gi = h.index("gcf_id")
            gni = h.index("genomes")
            for line in f:
                p = line.strip().split("\t")
                for g in p[gni].split(","):
                    if g.strip():
                        gcf_map[p[gi]].add(g.strip())

    niches = list(NICHE_PATTERNS.keys())
    print(f"Classes: {len(cls_map)}, GCFs: {len(gcf_map)}")

    for label, fmap, out in [("class", cls_map, "ecological_association_classes.tsv"),
                              ("gcf", gcf_map, "ecological_association_gcfs.tsv")]:
        print(f"\nTesting {label} x niche...")
        res = test_assoc(genomes, fmap, niches_by, niches)
        if res:
            _, fdr, _, _ = multipletests([r["p_value"] for r in res], method="fdr_bh")
            for r, v in zip(res, fdr):
                r["fdr"] = v
            res.sort(key=lambda r: r["p_value"])
            with open(RES_STATS / out, "w") as f:
                f.write(f"{label}\tniche\ta\tb\tc\td\todds_ratio\tp_value\tfdr\n")
                for r in res:
                    f.write(f"{r['feature']}\t{r['niche']}\t{r['a']}\t{r['b']}\t{r['c']}\t{r['d']}\t{r['odds_ratio']:.3f}\t{r['p_value']:.4e}\t{r['fdr']:.4e}\n")
            n_sig = sum(1 for r in res if r["fdr"] < 0.05)
            print(f"  Tests: {len(res)}, Sig (FDR<0.05): {n_sig}")


if __name__ == "__main__":
    main()
