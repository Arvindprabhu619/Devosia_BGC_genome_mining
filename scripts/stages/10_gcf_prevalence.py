#!/usr/bin/env python3
# ============================================================
# scripts/10_gcf_prevalence.py
# ============================================================
# Purpose: Compute GCF prevalence from BiG-SCAPE 2.0 SQLite DB.
# ============================================================
import sqlite3
import csv
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / "results" / "09_bigscape" / "raw_output" / "raw_output.db"
PARSED = ROOT / "results" / "09_bigscape" / "parsed"
GTDBTK = ROOT / "results" / "05_gtdbtk"

CORE_THRESHOLD = 90.0
RARE_THRESHOLD = 10.0


def path_to_acc(p):
    """Extract Devosia accession from a BiG-SCAPE input path."""
    name = Path(p).stem
    if "__" in name:
        return name.split("__")[0]
    return name.split("_region")[0]


def main():
    print("=" * 60)
    print("10_gcf_prevalence.py (SQLite)")
    print("=" * 60)
    PARSED.mkdir(parents=True, exist_ok=True)

    final_list = GTDBTK / "final_genome_list.txt"
    genomes = [g.strip() for g in final_list.read_text().splitlines() if g.strip()]
    n_genomes = len(genomes)
    genome_set = set(genomes)
    print(f"Genomes: {n_genomes}")

    devosia_a = set()
    tax_file = GTDBTK / "taxonomy_assignments.tsv"
    if tax_file.exists():
        with open(tax_file) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row.get("status") == "Devosia_A":
                    devosia_a.add(row["accession"])
    print(f"Devosia_A: {len(devosia_a)}")

    if not DB.exists():
        print(f"ERROR: {DB} not found")
        return
    con = sqlite3.connect(str(DB))
    cur = con.cursor()

    for cutoff in [0.3, 0.5, 0.7]:
        print(f"\n--- Cutoff c{cutoff} ---")
        cur.execute("""
            SELECT g.path, rf.family_id
            FROM bgc_record r
            JOIN gbk g ON r.gbk_id = g.id
            JOIN bgc_record_family rf ON rf.record_id = r.id
            JOIN family f ON rf.family_id = f.id
            WHERE r.record_type = 'region'
              AND f.cutoff = ?
              AND (g.path LIKE '%GCF_%' OR g.path LIKE '%GCA_%')
        """, (cutoff,))
        rows = cur.fetchall()
        print(f"  Total Devosia BGC assignments: {len(rows)}")

        gcf_genomes = defaultdict(set)
        for path, family_id in rows:
            acc = path_to_acc(path)
            if acc in genome_set:
                gcf_genomes[family_id].add(acc)

        print(f"  Unique GCFs: {len(gcf_genomes)}")

        out_rows = []
        for gcf_id, accs in gcf_genomes.items():
            n = len(accs)
            prev = 100.0 * n / n_genomes
            if prev >= CORE_THRESHOLD:
                cat = "core"
            elif prev < RARE_THRESHOLD:
                cat = "rare"
            else:
                cat = "intermediate"
            out_rows.append({
                "gcf_id": f"FAM_{gcf_id:05d}",
                "n_genomes": n,
                "prevalence": round(prev, 2),
                "category": cat,
                "genomes": ",".join(sorted(accs)),
            })

        out_rows.sort(key=lambda r: -r["prevalence"])
        out_file = PARSED / f"gcf_prevalence_c{cutoff}.tsv"
        with open(out_file, "w") as f:
            f.write("gcf_id\tn_genomes\tprevalence\tcategory\tgenomes\n")
            for r in out_rows:
                f.write(f"{r['gcf_id']}\t{r['n_genomes']}\t{r['prevalence']}\t{r['category']}\t{r['genomes']}\n")

        nc = sum(1 for r in out_rows if r["category"] == "core")
        ni = sum(1 for r in out_rows if r["category"] == "intermediate")
        nr = sum(1 for r in out_rows if r["category"] == "rare")
        print(f"  GCFs: {len(out_rows)} | core={nc}, intermediate={ni}, rare={nr}")

    # Devosia vs Devosia_A overlap at c0.7
    c07 = PARSED / "gcf_prevalence_c0.7.tsv"
    if c07.exists() and devosia_a:
        dev_only, da_only, shared = set(), set(), set()
        with open(c07) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                accs = set(row["genomes"].split(",")) if row["genomes"] else set()
                in_dev = bool(accs - devosia_a)
                in_da = bool(accs & devosia_a)
                if in_dev and in_da:
                    shared.add(row["gcf_id"])
                elif in_dev:
                    dev_only.add(row["gcf_id"])
                elif in_da:
                    da_only.add(row["gcf_id"])
        with open(PARSED / "gcf_overlap_devosia_A.tsv", "w") as f:
            f.write("category\tn_gcfs\n")
            f.write(f"Devosia_specific\t{len(dev_only)}\n")
            f.write(f"Devosia_A_specific\t{len(da_only)}\n")
            f.write(f"shared\t{len(shared)}\n")
        print(f"\n--- Devosia vs Devosia_A at c0.7 ---")
        print(f"  Devosia-specific:   {len(dev_only)}")
        print(f"  Devosia_A-specific: {len(da_only)}")
        print(f"  Shared:             {len(shared)}")

    con.close()
    print("\n=== [10] Done ===")


if __name__ == "__main__":
    main()
