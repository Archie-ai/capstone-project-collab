#!/usr/bin/env python3

"""
Canonical preprocessing for the MicrobiomeHD par_scheperjans dataset.

Project:
Microbiome Abundance Consensus Machine

Dataset:
Scheperjans et al. Parkinson's disease microbiome study
MicrobiomeHD identifier: par_scheperjans

Primary preprocessing decisions:
- H   -> Control
- PAR -> Case
- Preserve all 148 samples because raw QC passed.
- Remove taxa present in fewer than 10% of samples.
- With 148 samples, 10% prevalence requires presence in >=15 samples.
- No additional total-count filter is applied.
- Assign stable analysis IDs while preserving complete source taxonomy.

Outputs are intended to become the canonical inputs used by:
1. CLR/Wilcoxon
2. ALDEx2
3. ANCOM-BC2
4. DESeq2
"""

from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd


# ============================================================
# 1. CONFIGURATION
# ============================================================

DATASET = "par_scheperjans"

PREVALENCE_THRESHOLD = 0.10

GROUP_MAP = {
    "H": "Control",
    "PAR": "Case",
}

BASE_RAW_DIR = Path(
    "data/raw/par_scheperjans/"
    "par_scheperjans_results"
)

METADATA_FILE = (
    BASE_RAW_DIR /
    "par_scheperjans.metadata.txt"
)

COUNTS_FILE = (
    BASE_RAW_DIR /
    "RDP" /
    "par_scheperjans.otu_table.100.denovo.rdp_assigned"
)

OUTPUT_DIR = Path(
    "data/processed/par_scheperjans"
)

RESULT_QC_DIR = Path(
    "results/real/par_scheperjans/preprocessing"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

RESULT_QC_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. OUTPUT FILES
# ============================================================

PROCESSED_COUNTS_FILE = (
    OUTPUT_DIR /
    "par_scheperjans_counts.csv"
)

PROCESSED_METADATA_FILE = (
    OUTPUT_DIR /
    "par_scheperjans_metadata.csv"
)

TAXONOMY_MAPPING_FILE = (
    OUTPUT_DIR /
    "par_scheperjans_taxonomy_mapping.csv"
)

TAXON_QC_FILE = (
    OUTPUT_DIR /
    "par_scheperjans_taxon_qc.csv"
)

PREPROCESSING_SUMMARY_FILE = (
    OUTPUT_DIR /
    "par_scheperjans_preprocessing_summary.txt"
)

QC_SUMMARY_FILE = (
    OUTPUT_DIR /
    "par_scheperjans_qc_summary.txt"
)

MANIFEST_FILE = (
    OUTPUT_DIR /
    "par_scheperjans_manifest.json"
)


# ============================================================
# 3. HELPERS
# ============================================================

def sha256_file(path: Path) -> str:
    """Calculate SHA-256 for a file."""

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b""
        ):
            digest.update(block)

    return digest.hexdigest()


def parse_taxonomy(source_taxon: str) -> dict:
    """
    Parse the semicolon-delimited RDP taxonomy string.

    Example:
    k__Bacteria;p__Firmicutes;c__Clostridia;...
    """

    fields = {
        "Kingdom": "",
        "Phylum": "",
        "Class": "",
        "Order": "",
        "Family": "",
        "Genus": "",
        "Species": "",
        "Denovo_ID": "",
    }

    prefix_map = {
        "k__": "Kingdom",
        "p__": "Phylum",
        "c__": "Class",
        "o__": "Order",
        "f__": "Family",
        "g__": "Genus",
        "s__": "Species",
        "d__": "Denovo_ID",
    }

    for token in str(source_taxon).split(";"):

        token = token.strip()

        for prefix, field in prefix_map.items():

            if token.startswith(prefix):

                value = token[len(prefix):].strip()

                fields[field] = value

                break

    return fields


def write_lines(
    path: Path,
    lines
) -> None:

    path.write_text(
        "\n".join(
            str(line)
            for line in lines
        ),
        encoding="utf-8"
    )


# ============================================================
# 4. LOAD SOURCE DATA
# ============================================================

print("=" * 72)
print("PAR_SCHEPERJANS CANONICAL PREPROCESSING")
print("=" * 72)

print("\nLoading source data...")

metadata_raw = pd.read_csv(
    METADATA_FILE,
    sep="\t",
    dtype=str
)

counts_raw = pd.read_csv(
    COUNTS_FILE,
    sep="\t",
    index_col=0
)

metadata_raw = metadata_raw.rename(
    columns={
        metadata_raw.columns[0]: "SampleID"
    }
)

metadata_raw["SampleID"] = (
    metadata_raw["SampleID"]
    .astype(str)
    .str.strip()
)

metadata_raw["DiseaseState"] = (
    metadata_raw["DiseaseState"]
    .astype(str)
    .str.strip()
)

counts_raw.columns = (
    counts_raw.columns
    .astype(str)
    .str.strip()
)

counts_raw.index = (
    counts_raw.index
    .astype(str)
    .str.strip()
)


print(
    f"Raw metadata: "
    f"{metadata_raw.shape[0]} samples x "
    f"{metadata_raw.shape[1]} columns"
)

print(
    f"Raw counts: "
    f"{counts_raw.shape[0]} taxa x "
    f"{counts_raw.shape[1]} samples"
)


# ============================================================
# 5. VALIDATE SOURCE DATA
# ============================================================

if metadata_raw["SampleID"].duplicated().any():
    raise ValueError(
        "Duplicate metadata SampleIDs detected."
    )

if counts_raw.columns.duplicated().any():
    raise ValueError(
        "Duplicate sample columns detected in count table."
    )

if counts_raw.index.duplicated().any():
    raise ValueError(
        "Duplicate source taxon labels detected."
    )

observed_groups = set(
    metadata_raw["DiseaseState"]
)

expected_groups = set(
    GROUP_MAP
)

if observed_groups != expected_groups:
    raise ValueError(
        "Unexpected DiseaseState labels.\n"
        f"Observed: {sorted(observed_groups)}\n"
        f"Expected: {sorted(expected_groups)}"
    )

metadata_samples = set(
    metadata_raw["SampleID"]
)

count_samples = set(
    counts_raw.columns
)

if metadata_samples != count_samples:
    raise ValueError(
        "Metadata and count-table sample sets do not match."
    )


# ============================================================
# 6. CONVERT COUNTS TO VALID INTEGER MATRIX
# ============================================================

counts_numeric = counts_raw.apply(
    pd.to_numeric,
    errors="coerce"
)

if counts_numeric.isna().any().any():
    raise ValueError(
        "Non-numeric or missing count values detected."
    )

values = counts_numeric.to_numpy(
    dtype=float
)

if not np.isfinite(values).all():
    raise ValueError(
        "Non-finite count values detected."
    )

if (values < 0).any():
    raise ValueError(
        "Negative count values detected."
    )

if not np.allclose(
    values,
    np.round(values)
):
    raise ValueError(
        "Non-integer count values detected."
    )

counts_numeric = (
    counts_numeric
    .round()
    .astype("int64")
)


# ============================================================
# 7. STANDARDIZE SAMPLE ORDER AND GROUP LABELS
# ============================================================

metadata = metadata_raw.copy()

metadata["Group"] = (
    metadata["DiseaseState"]
    .map(GROUP_MAP)
)

if metadata["Group"].isna().any():
    raise ValueError(
        "One or more DiseaseState values could not be mapped."
    )


# Preserve original metadata order.

sample_order = (
    metadata["SampleID"]
    .tolist()
)

counts_numeric = counts_numeric[
    sample_order
]


if list(counts_numeric.columns) != sample_order:
    raise ValueError(
        "Count columns did not align to metadata order."
    )


# ============================================================
# 8. PRE-FILTER CHARACTERISTICS
# ============================================================

n_samples = counts_numeric.shape[1]

n_taxa_raw = counts_numeric.shape[0]

minimum_samples_required = int(
    math.ceil(
        PREVALENCE_THRESHOLD *
        n_samples
    )
)

present_counts_raw = (
    counts_numeric > 0
).sum(axis=1)

prevalence_raw = (
    present_counts_raw /
    n_samples
)

total_counts_raw = (
    counts_numeric.sum(axis=1)
)

zero_fraction_raw = float(
    (
        counts_numeric == 0
    ).sum().sum()
    /
    counts_numeric.size
)


# ============================================================
# 9. APPLY PRIMARY PREVALENCE FILTER
# ============================================================

retain_mask = (
    present_counts_raw >=
    minimum_samples_required
)

counts_filtered_source = (
    counts_numeric.loc[
        retain_mask
    ].copy()
)

n_taxa_retained = (
    counts_filtered_source.shape[0]
)

n_taxa_removed = (
    n_taxa_raw -
    n_taxa_retained
)

print("\nPrevalence filtering:")

print(
    f"Threshold: "
    f"{PREVALENCE_THRESHOLD:.0%}"
)

print(
    f"Samples: "
    f"{n_samples}"
)

print(
    f"Minimum samples required: "
    f"{minimum_samples_required}"
)

print(
    f"Raw taxa: "
    f"{n_taxa_raw}"
)

print(
    f"Retained taxa: "
    f"{n_taxa_retained}"
)

print(
    f"Removed taxa: "
    f"{n_taxa_removed}"
)


# ============================================================
# 10. BUILD STABLE ANALYSIS TAXON IDs
# ============================================================

analysis_ids = [
    f"OTU_{i:05d}"
    for i in range(
        1,
        n_taxa_retained + 1
    )
]

source_taxa_retained = (
    counts_filtered_source.index
    .astype(str)
    .tolist()
)

source_to_analysis = dict(
    zip(
        source_taxa_retained,
        analysis_ids
    )
)

counts_filtered = (
    counts_filtered_source.copy()
)

counts_filtered.index = [
    source_to_analysis[
        source_taxon
    ]
    for source_taxon
    in counts_filtered_source.index
]

counts_filtered.index.name = "Taxon"


# ============================================================
# 11. BUILD TAXONOMY MAPPING
# ============================================================

taxonomy_rows = []

for source_taxon in source_taxa_retained:

    parsed = parse_taxonomy(
        source_taxon
    )

    analysis_id = (
        source_to_analysis[
            source_taxon
        ]
    )

    row = {
        "Taxon": analysis_id,
        "Source_Taxon": source_taxon,
        "Denovo_ID": parsed[
            "Denovo_ID"
        ],
        "Kingdom": parsed[
            "Kingdom"
        ],
        "Phylum": parsed[
            "Phylum"
        ],
        "Class": parsed[
            "Class"
        ],
        "Order": parsed[
            "Order"
        ],
        "Family": parsed[
            "Family"
        ],
        "Genus": parsed[
            "Genus"
        ],
        "Species": parsed[
            "Species"
        ],
    }

    taxonomy_rows.append(row)


taxonomy_mapping = pd.DataFrame(
    taxonomy_rows
)


if taxonomy_mapping["Taxon"].duplicated().any():
    raise ValueError(
        "Duplicate analysis Taxon IDs generated."
    )

if taxonomy_mapping["Source_Taxon"].duplicated().any():
    raise ValueError(
        "Duplicate source taxa detected in taxonomy mapping."
    )

if taxonomy_mapping["Denovo_ID"].duplicated().any():
    raise ValueError(
        "Duplicate denovo identifiers detected after filtering."
    )

if taxonomy_mapping["Denovo_ID"].eq("").any():
    raise ValueError(
        "Missing denovo identifiers detected after filtering."
    )


# ============================================================
# 12. GROUP-SPECIFIC TAXON QC
# ============================================================

control_ids = (
    metadata.loc[
        metadata["Group"] ==
        "Control",
        "SampleID"
    ]
    .tolist()
)

case_ids = (
    metadata.loc[
        metadata["Group"] ==
        "Case",
        "SampleID"
    ]
    .tolist()
)

present_control = (
    counts_filtered[
        control_ids
    ] > 0
).sum(axis=1)

present_case = (
    counts_filtered[
        case_ids
    ] > 0
).sum(axis=1)

present_total = (
    counts_filtered > 0
).sum(axis=1)

prevalence_total = (
    present_total /
    n_samples
)

total_count_filtered = (
    counts_filtered.sum(axis=1)
)


taxon_qc = pd.DataFrame(
    {
        "Taxon":
            counts_filtered.index,

        "Total_Count":
            total_count_filtered.values,

        "Present_Samples":
            present_total.values,

        "Prevalence":
            prevalence_total.values,

        "Present_Control":
            present_control.values,

        "Control_Prevalence":
            (
                present_control /
                len(control_ids)
            ).values,

        "Present_Case":
            present_case.values,

        "Case_Prevalence":
            (
                present_case /
                len(case_ids)
            ).values,
    }
)

taxon_qc = taxon_qc.merge(
    taxonomy_mapping[
        [
            "Taxon",
            "Denovo_ID",
            "Phylum",
            "Family",
            "Genus",
        ]
    ],
    on="Taxon",
    how="left",
    validate="one_to_one"
)


# ============================================================
# 13. POST-FILTER QC
# ============================================================

post_library_sizes = (
    counts_filtered.sum(axis=0)
)

zero_library_samples_post = (
    post_library_sizes[
        post_library_sizes == 0
    ]
)

all_zero_taxa_post = int(
    (
        counts_filtered.sum(
            axis=1
        ) == 0
    ).sum()
)

zero_fraction_post = float(
    (
        counts_filtered == 0
    ).sum().sum()
    /
    counts_filtered.size
)

minimum_prevalence_post = int(
    (
        counts_filtered > 0
    ).sum(axis=1).min()
)

median_prevalence_post = float(
    (
        counts_filtered > 0
    ).sum(axis=1).median()
)

maximum_prevalence_post = int(
    (
        counts_filtered > 0
    ).sum(axis=1).max()
)


print("\nPost-filter QC:")

print(
    f"Filtered matrix: "
    f"{counts_filtered.shape[0]} taxa x "
    f"{counts_filtered.shape[1]} samples"
)

print(
    f"Control samples: "
    f"{len(control_ids)}"
)

print(
    f"Case samples: "
    f"{len(case_ids)}"
)

print(
    f"Zero fraction after filtering: "
    f"{zero_fraction_post:.6f}"
)

print(
    f"Zero percentage after filtering: "
    f"{zero_fraction_post * 100:.2f}%"
)

print(
    f"All-zero taxa after filtering: "
    f"{all_zero_taxa_post}"
)

print(
    f"Zero-library samples after filtering: "
    f"{len(zero_library_samples_post)}"
)

print(
    f"Minimum retained taxon prevalence: "
    f"{minimum_prevalence_post}/{n_samples}"
)


# ============================================================
# 14. BUILD CANONICAL METADATA
# ============================================================

metadata_processed = metadata[
    [
        "SampleID",
        "Group",
        "DiseaseState",
        "instrument_model",
        "scientific_name",
        "study_accession",
        "secondary_study_accession",
        "sample_accession",
        "run_accession",
        "fastq_ftp",
    ]
].copy()


# ============================================================
# 15. FINAL VALIDATION
# ============================================================

checks = {
    "148 samples retained":
        counts_filtered.shape[1] == 148,

    "74 Control samples":
        len(control_ids) == 74,

    "74 Case samples":
        len(case_ids) == 74,

    "Metadata/count sample order aligned":
        list(
            counts_filtered.columns
        )
        ==
        metadata_processed[
            "SampleID"
        ].tolist(),

    "No duplicate Taxon IDs":
        not counts_filtered.index.duplicated().any(),

    "No duplicate sample IDs":
        not counts_filtered.columns.duplicated().any(),

    "No negative counts":
        not (
            counts_filtered < 0
        ).any().any(),

    "All counts integer":
        np.allclose(
            counts_filtered.to_numpy(),
            np.round(
                counts_filtered.to_numpy()
            )
        ),

    "No all-zero retained taxa":
        all_zero_taxa_post == 0,

    "No zero-library retained samples":
        len(
            zero_library_samples_post
        ) == 0,

    "Every retained taxon meets prevalence threshold":
        minimum_prevalence_post
        >=
        minimum_samples_required,

    "Taxonomy mapping complete":
        len(taxonomy_mapping)
        ==
        len(counts_filtered),

    "Taxon QC complete":
        len(taxon_qc)
        ==
        len(counts_filtered),
}


print("\nCanonical dataset checks:")

for label, passed in checks.items():

    print(
        f"{label}: "
        f"{'PASS' if passed else 'FAIL'}"
    )


overall_pass = all(
    checks.values()
)

if not overall_pass:
    raise RuntimeError(
        "Canonical preprocessing failed one or more QC checks."
    )


# ============================================================
# 16. WRITE CANONICAL FILES
# ============================================================

counts_filtered.to_csv(
    PROCESSED_COUNTS_FILE
)

metadata_processed.to_csv(
    PROCESSED_METADATA_FILE,
    index=False
)

taxonomy_mapping.to_csv(
    TAXONOMY_MAPPING_FILE,
    index=False
)

taxon_qc.to_csv(
    TAXON_QC_FILE,
    index=False
)


# ============================================================
# 17. HASH CANONICAL OUTPUTS
# ============================================================

hashes = {
    "counts_sha256":
        sha256_file(
            PROCESSED_COUNTS_FILE
        ),

    "metadata_sha256":
        sha256_file(
            PROCESSED_METADATA_FILE
        ),

    "taxonomy_mapping_sha256":
        sha256_file(
            TAXONOMY_MAPPING_FILE
        ),

    "taxon_qc_sha256":
        sha256_file(
            TAXON_QC_FILE
        ),
}


# ============================================================
# 18. PREPROCESSING SUMMARY
# ============================================================

preprocessing_lines = [
    "PAR_SCHEPERJANS PREPROCESSING SUMMARY",
    "=" * 68,
    "",
    "DATASET",
    "MicrobiomeHD: par_scheperjans",
    "",
    "SOURCE DATA",
    f"Raw samples: {n_samples}",
    f"Raw taxa: {n_taxa_raw}",
    (
        f"Raw zero fraction: "
        f"{zero_fraction_raw:.6f}"
    ),
    (
        f"Raw zero percentage: "
        f"{zero_fraction_raw * 100:.2f}%"
    ),
    "",
    "GROUP MAPPING",
    "H -> Control",
    "PAR -> Case",
    (
        f"Control samples: "
        f"{len(control_ids)}"
    ),
    (
        f"Case samples: "
        f"{len(case_ids)}"
    ),
    "",
    "PRIMARY FILTER",
    (
        "Prevalence threshold: "
        f"{PREVALENCE_THRESHOLD:.0%}"
    ),
    (
        "Minimum samples required: "
        f"{minimum_samples_required}"
    ),
    (
        "Filtering criterion: "
        "taxon must have count > 0 in at least "
        f"{minimum_samples_required} of {n_samples} samples"
    ),
    "Additional total-count filter: none",
    "",
    "FILTER RESULT",
    (
        f"Taxa retained: "
        f"{n_taxa_retained}"
    ),
    (
        f"Taxa removed: "
        f"{n_taxa_removed}"
    ),
    (
        f"Percent retained: "
        f"{n_taxa_retained / n_taxa_raw * 100:.2f}%"
    ),
    "",
    "POST-FILTER DATASET",
    (
        f"Shape: "
        f"{counts_filtered.shape[0]} taxa x "
        f"{counts_filtered.shape[1]} samples"
    ),
    (
        f"Zero fraction: "
        f"{zero_fraction_post:.6f}"
    ),
    (
        f"Zero percentage: "
        f"{zero_fraction_post * 100:.2f}%"
    ),
    (
        f"Minimum taxon prevalence: "
        f"{minimum_prevalence_post}"
    ),
    (
        f"Median taxon prevalence: "
        f"{median_prevalence_post:.1f}"
    ),
    (
        f"Maximum taxon prevalence: "
        f"{maximum_prevalence_post}"
    ),
    "",
    "ANALYSIS IDENTIFIERS",
    (
        "Stable Taxon IDs assigned as "
        "OTU_00001 ... OTU_N."
    ),
    (
        "Original RDP taxonomy and d__denovo "
        "identifier preserved in taxonomy mapping."
    ),
]


write_lines(
    PREPROCESSING_SUMMARY_FILE,
    preprocessing_lines
)


# ============================================================
# 19. QC SUMMARY
# ============================================================

qc_lines = [
    "PAR_SCHEPERJANS CANONICAL DATASET QC",
    "=" * 68,
    "",
    (
        f"Processed shape: "
        f"{counts_filtered.shape[0]} taxa x "
        f"{counts_filtered.shape[1]} samples"
    ),
    (
        f"Control samples: "
        f"{len(control_ids)}"
    ),
    (
        f"Case samples: "
        f"{len(case_ids)}"
    ),
    "",
    "POST-FILTER LIBRARY SIZE",
    (
        f"Mean: "
        f"{post_library_sizes.mean():.2f}"
    ),
    (
        f"Median: "
        f"{post_library_sizes.median():.2f}"
    ),
    (
        f"Minimum: "
        f"{post_library_sizes.min()}"
    ),
    (
        f"Maximum: "
        f"{post_library_sizes.max()}"
    ),
    "",
    "SPARSITY",
    (
        f"Zero fraction: "
        f"{zero_fraction_post:.6f}"
    ),
    (
        f"Zero percentage: "
        f"{zero_fraction_post * 100:.2f}%"
    ),
    "",
    "CHECKS",
]

for label, passed in checks.items():

    qc_lines.append(
        f"{label}: "
        f"{'PASS' if passed else 'FAIL'}"
    )


qc_lines.extend(
    [
        "",
        "SHA-256",
        (
            f"counts: "
            f"{hashes['counts_sha256']}"
        ),
        (
            f"metadata: "
            f"{hashes['metadata_sha256']}"
        ),
        (
            "taxonomy_mapping: "
            f"{hashes['taxonomy_mapping_sha256']}"
        ),
        (
            f"taxon_qc: "
            f"{hashes['taxon_qc_sha256']}"
        ),
        "",
        "OVERALL RESULT: PASS",
    ]
)


write_lines(
    QC_SUMMARY_FILE,
    qc_lines
)


# ============================================================
# 20. JSON MANIFEST
# ============================================================

manifest = {
    "dataset": DATASET,
    "raw_metadata_file": str(
        METADATA_FILE
    ),
    "raw_counts_file": str(
        COUNTS_FILE
    ),
    "raw_samples": int(
        n_samples
    ),
    "raw_taxa": int(
        n_taxa_raw
    ),
    "prevalence_threshold": (
        PREVALENCE_THRESHOLD
    ),
    "minimum_samples_required": int(
        minimum_samples_required
    ),
    "retained_taxa": int(
        n_taxa_retained
    ),
    "removed_taxa": int(
        n_taxa_removed
    ),
    "control_samples": int(
        len(control_ids)
    ),
    "case_samples": int(
        len(case_ids)
    ),
    "raw_zero_fraction": (
        zero_fraction_raw
    ),
    "processed_zero_fraction": (
        zero_fraction_post
    ),
    "group_mapping": (
        GROUP_MAP
    ),
    "additional_total_count_filter": None,
    "hashes": hashes,
    "overall_qc": "PASS",
}


MANIFEST_FILE.write_text(
    json.dumps(
        manifest,
        indent=2
    ),
    encoding="utf-8"
)


# ============================================================
# 21. FINAL OUTPUT
# ============================================================

print("\n" + "=" * 72)
print("CANONICAL PREPROCESSING: PASS")
print("=" * 72)

print("\nCanonical files:")

print(
    PROCESSED_COUNTS_FILE
)

print(
    PROCESSED_METADATA_FILE
)

print(
    TAXONOMY_MAPPING_FILE
)

print(
    TAXON_QC_FILE
)

print(
    PREPROCESSING_SUMMARY_FILE
)

print(
    QC_SUMMARY_FILE
)

print(
    MANIFEST_FILE
)

print("\nSHA-256:")

for name, digest in hashes.items():
    print(
        f"{name}: {digest}"
    )
