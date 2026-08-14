from pathlib import Path

p = Path("09_build_final_prioritized_BGC_table.py")
s = p.read_text()

# Replace the coordinate resolver block, if present
start_marker = "# ============================================================\n# RESOLVE BGC COORDINATES"
end_marker = "# ============================================================\n# PUBLICATION TABLE"

start = s.find(start_marker)
end = s.find(end_marker)

if start != -1 and end != -1:
    s = s[:start] + s[end:]

# Replace direct coordinate references in publication table
s = s.replace(
    'rep["region_number"].astype(str)',
    'rep["region_number"].astype(str)'
)

s = s.replace(
    'rep["start"].astype(str)',
    'rep["start_x"].astype(str)'
)

s = s.replace(
    'rep["end"].astype(str)',
    'rep["end_x"].astype(str)'
)

p.write_text(s)

print("Coordinate references corrected:")
print("  region -> region_number")
print("  start  -> start_x")
print("  end    -> end_x")
