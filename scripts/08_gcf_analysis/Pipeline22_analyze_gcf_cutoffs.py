import re
path = "analyze_gcf_cutoffs.py"
with open(path) as f:
    content = f.read()

old = '''    all_records = pd.concat(dfs, ignore_index=True)
    all_records[["contig", "region_number"]] = all_records["GBK"].apply(
        lambda x: pd.Series(parse_gbk(x))
    )'''

new = '''    all_records = pd.concat(dfs, ignore_index=True)
    n_before_filter = len(all_records)
    all_records = all_records[~all_records["GBK"].str.match(r"^BGC\\d{7}$")].copy()
    n_mibig_removed = n_before_filter - len(all_records)
    if n_mibig_removed:
        print(f"  (cutoff {cutoff}: removed {n_mibig_removed} MIBiG reference records)")
    all_records[["contig", "region_number"]] = all_records["GBK"].apply(
        lambda x: pd.Series(parse_gbk(x))
    )'''

assert old in content, "pattern not found -- check script hasn't changed"
content = content.replace(old, new)
with open(path, "w") as f:
    f.write(content)
print("patched OK")
