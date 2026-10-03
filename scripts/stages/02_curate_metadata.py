#!/usr/bin/env python3
# ============================================================
# scripts/02_curate_metadata.py
# ============================================================
# Purpose: Filter NCBI metadata to keep only cultured, non-fragmented
#          Devosia genomes.
#
# Inputs:   data/metadata/devosia_refseq.jsonl
#           data/metadata/devosia_genbank.jsonl
# Outputs:
#   data/metadata/curated_accessions.txt
#   data/metadata/curation_report.tsv
#   data/metadata/curation_summary.txt
#
# Removes:  MAGs, uncultured/Candidatus, >200 scaffolds
# Runtime:  2 min
# ============================================================
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata"
MAX_SCAFFOLDS = 200

MAG_PATTERNS = [
    r"\bMAG\b", r"\bMetaBAT\b", r"\bCONCOCT\b", r"\bMaxBin\b", r"\bBinChicken\b",
    r"\bbin_\b", r"\bbin\.\b", r"\bbinner\b", r"\brefined\b", r"^bin_",
    r"_bin_", r"\bmetagenom", r"\bmetagenome-assembled\b", r"\bUBA\d+",
    r"\bbin\b.*\bUBC\b",
]
MAG_RE = re.compile("|".join(MAG_PATTERNS), re.IGNORECASE)
UNCULT_RE = re.compile(r"\bCandidatus\b|\buncultured\b|\bunclassified\b", re.IGNORECASE)


def load_jsonl(p):
    out = []
    if not p.exists():
        return out
    with open(p) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return out


def extract(rec):
    acc = rec.get("accession", "")
    org = rec.get("organism", {}) or {}
    org_name = org.get("organism_name", "")
    infraspecific = org.get("infraspecific_names", {}) or {}
    strain = infraspecific.get("strain", "") or infraspecific.get("isolate", "")
    ai = rec.get("assembly_info", {}) or {}
    aname = ai.get("assembly_name", "")
    stats = rec.get("assembly_stats", {}) or {}
    n_scaf = stats.get("number_of_scaffolds", None)
    bs = rec.get("biosample", {}) or {}
    attrs = bs.get("attributes", []) or []
    bs_dict = {a.get("name", "").lower(): a.get("value", "") for a in attrs if a.get("name")}
    return {
        "accession": acc,
        "organism_name": org_name,
        "strain": strain,
        "assembly_name": aname,
        "n_scaffolds": n_scaf,
        "biosample": bs_dict,
    }


def main():
    print("=" * 60)
    print("02_curate_metadata.py")
    print("=" * 60)

    refseq = load_jsonl(META / "devosia_refseq.jsonl")
    genbank = load_jsonl(META / "devosia_genbank.jsonl")
    all_recs = refseq + genbank
    print(f"Loaded: RefSeq={len(refseq)}, GenBank={len(genbank)}, Total={len(all_recs)}")

    rows = [extract(r) for r in all_recs]

    # Deduplicate GenBank/RefSeq pairs
    def key(r):
        return (r["organism_name"].lower(), r["strain"].lower(), r["n_scaffolds"])

    deduped = {}
    for r in rows:
        k = key(r)
        is_refseq = r["accession"].startswith("GCF_")
        if k not in deduped:
            deduped[k] = r
        else:
            if is_refseq and not deduped[k]["accession"].startswith("GCF_"):
                deduped[k] = r
    deduped = list(deduped.values())
    print(f"After dedup: {len(rows)} -> {len(deduped)}")

    kept = []
    report = ["accession\torganism\tstrain\tassembly_name\tn_scaffolds\tdecision\treason"]

    for r in deduped:
        acc = r["accession"]
        text = f"{r['organism_name']} {r['strain']} {r['assembly_name']}"

        if MAG_RE.search(text):
            report.append(f"{acc}\t{r['organism_name']}\t{r['strain']}\t{r['assembly_name']}\t{r['n_scaffolds']}\tEXCLUDE_MAG\tkeyword")
            continue

        bs_text = " ".join(r["biosample"].values())
        if MAG_RE.search(bs_text) or "metagenome" in bs_text.lower():
            report.append(f"{acc}\t{r['organism_name']}\t{r['strain']}\t{r['assembly_name']}\t{r['n_scaffolds']}\tEXCLUDE_MAG\tbiosample")
            continue

        if UNCULT_RE.search(text):
            report.append(f"{acc}\t{r['organism_name']}\t{r['strain']}\t{r['assembly_name']}\t{r['n_scaffolds']}\tEXCLUDE_UNCULTURED\tname")
            continue

        if r["n_scaffolds"] is not None and r["n_scaffolds"] > MAX_SCAFFOLDS:
            report.append(f"{acc}\t{r['organism_name']}\t{r['strain']}\t{r['assembly_name']}\t{r['n_scaffolds']}\tEXCLUDE_FRAGMENTED\t>200 scaffolds")
            continue

        report.append(f"{acc}\t{r['organism_name']}\t{r['strain']}\t{r['assembly_name']}\t{r['n_scaffolds']}\tKEEP\tpassed")
        kept.append(r)

    (META / "curation_report.tsv").write_text("\n".join(report))
    (META / "curated_accessions.txt").write_text("\n".join(r["accession"] for r in kept) + "\n")

    n_mag = sum(1 for l in report if "EXCLUDE_MAG" in l)
    n_unc = sum(1 for l in report if "EXCLUDE_UNCULTURED" in l)
    n_frag = sum(1 for l in report if "EXCLUDE_FRAGMENTED" in l)

    summary = f"""=== Curation Summary ===
Raw records:                {len(all_recs)}
After dedup:                {len(deduped)}
Excluded (MAG):             {n_mag}
Excluded (uncultured):      {n_unc}
Excluded (fragmented):      {n_frag}
----------------------------------
FINAL curated assemblies:   {len(kept)}
"""
    print(summary)
    (META / "curation_summary.txt").write_text(summary)


if __name__ == "__main__":
    main()
