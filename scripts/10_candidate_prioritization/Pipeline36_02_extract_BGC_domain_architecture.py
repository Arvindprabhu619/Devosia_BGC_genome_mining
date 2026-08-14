#!/usr/bin/env python3

import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO


MASTER = "BGC_candidate_master_c0.7.tsv"
ANTISMASH_ROOT = Path("antismash_output")

OUT_BGC = "BGC_domain_architecture_c0.7.tsv"
OUT_DOMAIN = "BGC_domain_features_c0.7.tsv"


# ============================================================
# HELPERS
# ============================================================

def parse_region_filename(path):
    """
    Example:
      NZ_JNNO01000034.1.region001.gbk

    Returns:
      contig = NZ_JNNO01000034.1
      region = 1
    """
    m = re.match(
        r"^(.+)\.region(\d+)\.gbk$",
        path.name
    )

    if not m:
        return None, None

    return m.group(1), int(m.group(2))


def qvalues(feature, key):
    vals = feature.qualifiers.get(key, [])

    if isinstance(vals, str):
        return [vals]

    return [str(x) for x in vals]


def unique_keep_order(values):
    out = []
    seen = set()

    for x in values:
        x = str(x).strip()

        if not x:
            continue

        if x not in seen:
            seen.add(x)
            out.append(x)

    return out


# ============================================================
# LOAD MASTER
# ============================================================

master = pd.read_csv(
    MASTER,
    sep="\t"
)

print("=" * 100)
print("DEVOSIA BGC DOMAIN ARCHITECTURE EXTRACTION")
print("=" * 100)

print("Master BGCs:", len(master))


# ============================================================
# EXACT REGION-LEVEL ANTI-SMASH FILES
# ============================================================

region_files = sorted(
    ANTISMASH_ROOT.rglob("*.region*.gbk")
)

print("Region GBKs:", len(region_files))

if len(region_files) == 0:
    raise SystemExit(
        "ERROR: no region-level *.region*.gbk files found"
    )


# ============================================================
# PARSE
# ============================================================

bgc_rows = []
domain_rows = []

seen_bgc = set()

for i, path in enumerate(region_files, 1):

    contig, region_number = parse_region_filename(path)

    if contig is None:
        print("WARNING: could not parse:", path)
        continue

    # Genome accession is the antiSMASH parent directory
    accession = path.parent.name

    bgc_key = (
        str(accession),
        str(contig),
        int(region_number)
    )

    if bgc_key in seen_bgc:
        print("WARNING duplicate:", path)
        continue

    seen_bgc.add(bgc_key)

    try:
        records = list(
            SeqIO.parse(
                str(path),
                "genbank"
            )
        )
    except Exception as e:
        print(
            "WARNING parse failure:",
            path,
            "|",
            repr(e)
        )
        continue

    if len(records) != 1:
        print(
            "WARNING unexpected record count:",
            len(records),
            path
        )
        continue

    record = records[0]

    # --------------------------------------------------------
    # Region-level information
    # --------------------------------------------------------

    region_start = None
    region_end = None

    region_products = []
    protocluster_products = []
    candidate_products = []
    candidate_numbers = []

    for feature in record.features:

        if feature.type == "region":

            region_start = int(
                feature.location.start
            ) + 1

            region_end = int(
                feature.location.end
            )

            region_products.extend(
                qvalues(feature, "product")
            )

            candidate_numbers.extend(
                qvalues(
                    feature,
                    "candidate_cluster_numbers"
                )
            )

        elif feature.type == "protocluster":

            protocluster_products.extend(
                qvalues(feature, "product")
            )

        elif feature.type in [
            "cand_cluster",
            "candidate_cluster"
        ]:

            candidate_products.extend(
                qvalues(feature, "product")
            )

            candidate_numbers.extend(
                qvalues(
                    feature,
                    "candidate_cluster_numbers"
                )
            )

    # --------------------------------------------------------
    # CDS / domain features
    # --------------------------------------------------------

    cds_count = 0

    cds_with_gene_functions = 0

    asdomain_count = 0
    pfam_count = 0

    asdomain_annotations = []
    pfam_annotations = []

    gene_functions = []

    domain_order = []

    biosynthetic_loci = set()

    for feature in record.features:

        q = feature.qualifiers

        if feature.type == "CDS":

            cds_count += 1

            gf = qvalues(
                feature,
                "gene_functions"
            )

            if gf:
                cds_with_gene_functions += 1

            gene_functions.extend(gf)

        elif feature.type == "aSDomain":

            asdomain_count += 1

            ads = qvalues(
                feature,
                "aSDomain"
            )

            locus = qvalues(
                feature,
                "locus_tag"
            )

            locus = (
                locus[0]
                if locus
                else ""
            )

            for annotation in ads:

                asdomain_annotations.append(
                    annotation
                )

                biosynthetic_loci.add(
                    locus
                )

                domain_order.append(
                    (
                        int(
                            feature.location.start
                        ) + 1,
                        int(
                            feature.location.end
                        ),
                        annotation,
                        locus
                    )
                )

                domain_rows.append({
                    "accession": accession,
                    "contig": contig,
                    "region_number": region_number,
                    "feature_type": "aSDomain",
                    "locus_tag": locus,
                    "annotation": annotation,
                    "start": int(
                        feature.location.start
                    ) + 1,
                    "end": int(
                        feature.location.end
                    ),
                    "source_file": str(path)
                })

        elif feature.type == "PFAM_domain":

            pfam_count += 1

            # Retain ALL useful PFAM qualifiers without
            # inventing interpretation.
            vals = []

            for key in [
                "pfam",
                "PFAM",
                "description",
                "label",
                "name"
            ]:
                vals.extend(
                    qvalues(feature, key)
                )

            if not vals:
                vals = ["PFAM_domain"]

            for annotation in vals:

                pfam_annotations.append(
                    annotation
                )

                domain_rows.append({
                    "accession": accession,
                    "contig": contig,
                    "region_number": region_number,
                    "feature_type": "PFAM_domain",
                    "locus_tag": (
                        qvalues(
                            feature,
                            "locus_tag"
                        ) or [""]
                    )[0],
                    "annotation": annotation,
                    "start": int(
                        feature.location.start
                    ) + 1,
                    "end": int(
                        feature.location.end
                    ),
                    "source_file": str(path)
                })

    # --------------------------------------------------------
    # Ordered aSDomain architecture
    # --------------------------------------------------------

    domain_order.sort(
        key=lambda x: (x[0], x[1])
    )

    ordered_domains = [
        x[2]
        for x in domain_order
    ]

    ordered_loci = [
        x[3]
        for x in domain_order
    ]

    # --------------------------------------------------------
    # Clean annotations
    # --------------------------------------------------------

    asdomain_annotations = unique_keep_order(
        asdomain_annotations
    )

    pfam_annotations = unique_keep_order(
        pfam_annotations
    )

    gene_functions = unique_keep_order(
        gene_functions
    )

    products = unique_keep_order(
        region_products
        + protocluster_products
        + candidate_products
    )

    candidate_numbers = unique_keep_order(
        candidate_numbers
    )

    # --------------------------------------------------------
    # Descriptive architecture flags
    #
    # IMPORTANT:
    # These are evidence flags, NOT biological scores.
    # --------------------------------------------------------

    architecture_flags = []

    product_text = " ".join(
        products
    ).lower()

    gf_text = " ".join(
        gene_functions
    ).lower()

    # Hybrid product assignment
    if ";" in product_text:
        architecture_flags.append(
            "multi_product_region"
        )

    # Multiple candidate clusters
    if len(candidate_numbers) > 1:
        architecture_flags.append(
            "multiple_candidate_clusters"
        )

    # Multiple aSDomains
    if len(asdomain_annotations) >= 4:
        architecture_flags.append(
            "multiple_asdomains"
        )

    if len(asdomain_annotations) >= 8:
        architecture_flags.append(
            "high_asdomain_count"
        )

    # Rule-based biosynthetic annotation
    if gene_functions:
        architecture_flags.append(
            "rule_based_biosynthetic_annotation"
        )

    # Search actual gene-function text only.
    #
    # We deliberately do NOT infer PKS/NRPS from arbitrary
    # TIGR IDs because the IDs alone are not sufficient.
    #

    if any(
        x in gf_text
        for x in [
            "polyketide",
            "pks",
            "ketosynthase",
            "acyltransferase"
        ]
    ):
        architecture_flags.append(
            "PKS_function_annotation"
        )

    if any(
        x in gf_text
        for x in [
            "nonribosomal",
            "nrps",
            "adenylation",
            "condensation"
        ]
    ):
        architecture_flags.append(
            "NRPS_function_annotation"
        )

    if any(
        x in gf_text
        for x in [
            "ripp",
            "ribosomally",
            "lanthipeptide",
            "lassopeptide"
        ]
    ):
        architecture_flags.append(
            "RiPP_function_annotation"
        )

    # --------------------------------------------------------
    # BGC-level row
    # --------------------------------------------------------

    bgc_rows.append({

        "accession": accession,

        "contig": contig,

        "region_number": region_number,

        "gbk_region_start": region_start,

        "gbk_region_end": region_end,

        "gbk_region_length": (
            region_end - region_start + 1
            if region_start is not None
            and region_end is not None
            else None
        ),

        "cds_count": cds_count,

        "biosynthetic_cds_count": max(
            len(biosynthetic_loci),
            cds_with_gene_functions
        ),

        "aSDomain_count": asdomain_count,

        "unique_aSDomain_count": len(
            asdomain_annotations
        ),

        "PFAM_feature_count": pfam_count,

        "unique_PFAM_annotation_count": len(
            pfam_annotations
        ),

        "aSDomain_annotations": ";".join(
            asdomain_annotations
        ),

        "ordered_aSDomain_architecture": " > ".join(
            ordered_domains
        ),

        "ordered_aSDomain_loci": ";".join(
            ordered_loci
        ),

        "PFAM_annotations": ";".join(
            pfam_annotations
        ),

        "gene_functions": ";".join(
            gene_functions
        ),

        "antiSMASH_products_from_GBK": ";".join(
            products
        ),

        "candidate_cluster_numbers_GBK": ";".join(
            candidate_numbers
        ),

        "architecture_flags": "; ".join(
            architecture_flags
        ),

        "source_file": str(path)
    })

    if i % 100 == 0:
        print(
            f"Processed {i}/{len(region_files)}"
        )


# ============================================================
# DATA FRAMES
# ============================================================

arch = pd.DataFrame(
    bgc_rows
)

domains = pd.DataFrame(
    domain_rows
)

# ============================================================
# MERGE WITH 582-BGC MASTER
# ============================================================

master2 = master.merge(
    arch,
    on=[
        "accession",
        "contig",
        "region_number"
    ],
    how="left",
    validate="one_to_one"
)

# ============================================================
# VALIDATION
# ============================================================

print("\n" + "=" * 100)
print("ARCHITECTURE VALIDATION")
print("=" * 100)

print(
    "Region-level files found:",
    len(region_files)
)

print(
    "Unique region files parsed:",
    len(arch)
)

print(
    "Master BGCs:",
    len(master)
)

print(
    "Merged BGCs:",
    len(master2)
)

print(
    "Architecture recovered:",
    master2["cds_count"].notna().sum()
)

print(
    "Architecture missing:",
    master2["cds_count"].isna().sum()
)

print(
    "aSDomain annotations present:",
    master2["aSDomain_count"].notna().sum()
)

print(
    "PFAM feature annotations present:",
    master2["PFAM_feature_count"].notna().sum()
)

# ============================================================
# IMPORTANT EXPECTED RESULT
# ============================================================

if len(arch) != 582:
    print(
        "\nWARNING: expected exactly 582 unique region records."
    )

if (
    master2["cds_count"].notna().sum()
    != 582
):
    print(
        "\nWARNING: architecture did not map to all 582 BGCs."
    )

# ============================================================
# SUMMARY
# ============================================================

print("\nCDS count:")
print(
    master2["cds_count"].describe().to_string()
)

print("\naSDomain count:")
print(
    master2["aSDomain_count"].describe().to_string()
)

print("\nPFAM feature count:")
print(
    master2["PFAM_feature_count"].describe().to_string()
)

print("\nArchitecture flags:")
print(
    master2["architecture_flags"]
    .fillna("")
    .replace("", "none")
    .value_counts()
    .head(30)
    .to_string()
)

print("\nTop 30 by aSDomain count:")
cols = [
    "accession",
    "contig",
    "region_number",
    "products",
    "novelty_category",
    "Family",
    "GCF_prevalence_pct",
    "GCF_category",
    "cds_count",
    "biosynthetic_cds_count",
    "aSDomain_count",
    "unique_aSDomain_count",
    "architecture_flags",
    "ordered_aSDomain_architecture"
]

print(
    master2.sort_values(
        [
            "aSDomain_count",
            "unique_aSDomain_count"
        ],
        ascending=False,
        na_position="last"
    )[cols]
    .head(30)
    .to_string(index=False)
)

# ============================================================
# SAVE
# ============================================================

master2.to_csv(
    OUT_BGC,
    sep="\t",
    index=False
)

domains.to_csv(
    OUT_DOMAIN,
    sep="\t",
    index=False
)

print("\n" + "=" * 100)
print("FILES GENERATED")
print("=" * 100)

print(OUT_BGC)
print(OUT_DOMAIN)
