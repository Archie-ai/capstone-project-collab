#!/usr/bin/env python3

"""
Raw-data audit for the MicrobiomeHD par_scheperjans dataset.

Project:
Microbiome Abundance Consensus Machine

Purpose:
- Validate metadata and OTU-table alignment.
- Confirm Parkinson's and healthy-control group counts.
- Validate raw count integrity.
- Quantify sequencing-depth characteristics.
- Quantify sparsity and taxon prevalence.
- Compare candidate prevalence-filter thresholds before
  committing to a canonical processed dataset.

No filtering or statistical testing is performed here.
"""

from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# 1. INPUTS
# ============================================================

BASE_DIR = Path(
    "data/raw/par_scheperjans/"
    "par_scheperjans_results"
)

METADATA_FILE = (
    BASE_DIR /
    "par_scheperjans.metadata.txt"
)

COUNTS_FILE = (
    BASE_DIR /
    "RDP/"
    "par_scheperjans.otu_table.100.denovo.rdp_assigned"
)

OUTPUT_DIR = Path(
    "results/real/par_scheperjans/raw_qc"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. LOAD DATA
# ============================================================

metadata = pd.read_csv(
    METADATA_FILE,
    sep="\t",
    dtype=str
)

counts = pd.read_csv(
    COUNTS_FILE,
    sep="\t",
    index_col=0
)


# ============================================================
# 3. STANDARDIZE METADATA SAMPLE-ID COLUMN
# ============================================================

metadata = metadata.rename(
    columns={
        metadata.columns[0]: "SampleID"
    }
)

metadata["SampleID"] = (
    metadata["SampleID"]
    .astype(str)
    .str.strip()
)

metadata["DiseaseState"] = (
    metadata["DiseaseState"]
    .astype(str)
    .str.strip()
)


# ============================================================
# 4. BASIC DIMENSIONS
# ============================================================

print("=" * 72)
print("PAR_SCHEPERJANS RAW DATA AUDIT")
print("=" * 72)

print("\nInput files:")
print(f"Metadata: {METADATA_FILE}")
print(f"Counts:   {COUNTS_FILE}")

print("\nRaw dimensions:")
print(
    f"Metadata: {metadata.shape[0]} samples x "
    f"{metadata.shape[1]} columns"
)

print(
    f"OTU table: {counts.shape[0]} taxa x "
    f"{counts.shape[1]} samples"
)


# ============================================================
# 5. METADATA QC
# ============================================================

duplicate_metadata_ids = (
    metadata["SampleID"]
    .duplicated()
    .sum()
)

missing_sample_ids = (
    metadata["SampleID"]
    .isna()
    .sum()
)

missing_disease_state = (
    metadata["DiseaseState"]
    .isna()
    .sum()
)

print("\nMetadata QC:")

print(
    f"Duplicate SampleIDs: "
    f"{duplicate_metadata_ids}"
)

print(
    f"Missing SampleIDs: "
    f"{missing_sample_ids}"
)

print(
    f"Missing DiseaseState values: "
    f"{missing_disease_state}"
)

print("\nDiseaseState counts:")

group_counts = (
    metadata["DiseaseState"]
    .value_counts(dropna=False)
)

print(
    group_counts.to_string()
)


# ============================================================
# 6. EXPECTED GROUP LABELS
# ============================================================

observed_groups = set(
    metadata["DiseaseState"].dropna()
)

expected_groups = {
    "H",
    "PAR",
}

unexpected_groups = (
    observed_groups -
    expected_groups
)

missing_expected_groups = (
    expected_groups -
    observed_groups
)

print("\nDiseaseState label validation:")

print(
    f"Observed groups: "
    f"{sorted(observed_groups)}"
)

print(
    f"Unexpected groups: "
    f"{sorted(unexpected_groups)}"
)

print(
    f"Missing expected groups: "
    f"{sorted(missing_expected_groups)}"
)


# ============================================================
# 7. SAMPLE-ID ALIGNMENT
# ============================================================

metadata_ids = set(
    metadata["SampleID"]
)

count_ids = set(
    counts.columns.astype(str)
)

metadata_not_counts = sorted(
    metadata_ids -
    count_ids
)

counts_not_metadata = sorted(
    count_ids -
    metadata_ids
)

print("\nSample alignment:")

print(
    f"Metadata sample IDs: "
    f"{len(metadata_ids)}"
)

print(
    f"Count-table sample IDs: "
    f"{len(count_ids)}"
)

print(
    f"Metadata IDs absent from counts: "
    f"{len(metadata_not_counts)}"
)

print(
    f"Count IDs absent from metadata: "
    f"{len(counts_not_metadata)}"
)

if metadata_not_counts:
    print(
        "\nMetadata-only IDs:"
    )
    print(
        "\n".join(metadata_not_counts)
    )

if counts_not_metadata:
    print(
        "\nCounts-only IDs:"
    )
    print(
        "\n".join(counts_not_metadata)
    )

sample_sets_match = (
    metadata_ids ==
    count_ids
)

print(
    f"\nSample sets match exactly: "
    f"{sample_sets_match}"
)


# ============================================================
# 8. COUNT-TABLE NUMERIC VALIDATION
# ============================================================

counts_numeric = counts.apply(
    pd.to_numeric,
    errors="coerce"
)

non_numeric_cells = int(
    counts_numeric
    .isna()
    .sum()
    .sum()
)

finite_mask = np.isfinite(
    counts_numeric.to_numpy(
        dtype=float
    )
)

non_finite_cells = int(
    (~finite_mask).sum()
)

negative_cells = int(
    (
        counts_numeric < 0
    )
    .sum()
    .sum()
)

values = counts_numeric.to_numpy(
    dtype=float
)

integer_mask = (
    np.isnan(values) |
    np.isclose(
        values,
        np.round(values)
    )
)

non_integer_cells = int(
    (~integer_mask).sum()
)

print("\nCount integrity:")

print(
    f"Non-numeric cells: "
    f"{non_numeric_cells}"
)

print(
    f"Non-finite cells: "
    f"{non_finite_cells}"
)

print(
    f"Negative cells: "
    f"{negative_cells}"
)

print(
    f"Non-integer cells: "
    f"{non_integer_cells}"
)


# ============================================================
# 9. STOP IF COUNTS ARE INVALID
# ============================================================

if non_numeric_cells > 0:
    raise ValueError(
        "Non-numeric values detected in OTU table."
    )

if non_finite_cells > 0:
    raise ValueError(
        "Non-finite values detected in OTU table."
    )

if negative_cells > 0:
    raise ValueError(
        "Negative values detected in OTU table."
    )

if non_integer_cells > 0:
    raise ValueError(
        "Non-integer values detected in OTU table."
    )


counts_numeric = (
    counts_numeric
    .round()
    .astype("int64")
)


# ============================================================
# 10. LIBRARY-SIZE QC
# ============================================================

library_sizes = (
    counts_numeric.sum(axis=0)
)

zero_library_samples = (
    library_sizes[
        library_sizes == 0
    ]
)

print("\nLibrary-size QC:")

print(
    f"Mean library size: "
    f"{library_sizes.mean():.2f}"
)

print(
    f"Median library size: "
    f"{library_sizes.median():.2f}"
)

print(
    f"Minimum library size: "
    f"{library_sizes.min()}"
)

print(
    f"Maximum library size: "
    f"{library_sizes.max()}"
)

print(
    f"Zero-library samples: "
    f"{len(zero_library_samples)}"
)


# ============================================================
# 11. SPARSITY
# ============================================================

total_cells = (
    counts_numeric.shape[0] *
    counts_numeric.shape[1]
)

zero_cells = int(
    (
        counts_numeric == 0
    )
    .sum()
    .sum()
)

zero_fraction = (
    zero_cells /
    total_cells
)

print("\nSparsity:")

print(
    f"Total taxon-sample cells: "
    f"{total_cells}"
)

print(
    f"Zero cells: "
    f"{zero_cells}"
)

print(
    f"Zero fraction: "
    f"{zero_fraction:.6f}"
)

print(
    f"Zero percentage: "
    f"{zero_fraction * 100:.2f}%"
)


# ============================================================
# 12. TAXON PREVALENCE
# ============================================================

present_samples = (
    counts_numeric > 0
).sum(axis=1)

prevalence_fraction = (
    present_samples /
    counts_numeric.shape[1]
)

all_zero_taxa = (
    present_samples == 0
)

print("\nTaxon prevalence:")

print(
    f"All-zero taxa: "
    f"{int(all_zero_taxa.sum())}"
)

print(
    f"Minimum prevalence count: "
    f"{int(present_samples.min())}"
)

print(
    f"Median prevalence count: "
    f"{present_samples.median():.1f}"
)

print(
    f"Maximum prevalence count: "
    f"{int(present_samples.max())}"
)

print(
    f"Median prevalence fraction: "
    f"{prevalence_fraction.median():.4f}"
)


# ============================================================
# 13. PREVALENCE-THRESHOLD COMPARISON
# ============================================================

thresholds = [
    0.01,
    0.05,
    0.10,
    0.15,
    0.20,
    0.25,
]

threshold_rows = []

n_samples = counts_numeric.shape[1]

for threshold in thresholds:

    minimum_samples = int(
        np.ceil(
            threshold * n_samples
        )
    )

    retained = (
        present_samples >=
        minimum_samples
    )

    retained_count = int(
        retained.sum()
    )

    removed_count = int(
        len(retained) -
        retained_count
    )

    threshold_rows.append(
        {
            "Prevalence_Threshold": threshold,
            "Minimum_Samples_Required": (
                minimum_samples
            ),
            "Taxa_Retained": retained_count,
            "Taxa_Removed": removed_count,
            "Percent_Taxa_Retained": (
                retained_count /
                len(retained) *
                100
            ),
        }
    )

threshold_df = pd.DataFrame(
    threshold_rows
)

print(
    "\nCandidate prevalence thresholds:"
)

print(
    threshold_df.to_string(
        index=False,
        formatters={
            "Prevalence_Threshold":
                lambda x: f"{x:.0%}",
            "Percent_Taxa_Retained":
                lambda x: f"{x:.2f}%",
        }
    )
)


# ============================================================
# 14. GROUP-SPECIFIC PREVALENCE
# ============================================================

metadata_indexed = (
    metadata
    .set_index("SampleID")
)

control_ids = (
    metadata_indexed.index[
        metadata_indexed[
            "DiseaseState"
        ] == "H"
    ]
    .tolist()
)

case_ids = (
    metadata_indexed.index[
        metadata_indexed[
            "DiseaseState"
        ] == "PAR"
    ]
    .tolist()
)

control_ids = [
    sample
    for sample in control_ids
    if sample in counts_numeric.columns
]

case_ids = [
    sample
    for sample in case_ids
    if sample in counts_numeric.columns
]

control_prevalence = (
    counts_numeric[
        control_ids
    ] > 0
).sum(axis=1)

case_prevalence = (
    counts_numeric[
        case_ids
    ] > 0
).sum(axis=1)

print("\nGroup-specific sample counts:")

print(
    f"Healthy/control samples: "
    f"{len(control_ids)}"
)

print(
    f"Parkinson's/case samples: "
    f"{len(case_ids)}"
)


# ============================================================
# 15. TAXON-LEVEL QC TABLE
# ============================================================

taxon_qc = pd.DataFrame(
    {
        "Source_Taxon": (
            counts_numeric.index
            .astype(str)
        ),
        "Total_Count": (
            counts_numeric.sum(axis=1)
            .to_numpy()
        ),
        "Present_Samples": (
            present_samples.to_numpy()
        ),
        "Prevalence": (
            prevalence_fraction.to_numpy()
        ),
        "Present_Healthy": (
            control_prevalence.to_numpy()
        ),
        "Present_PAR": (
            case_prevalence.to_numpy()
        ),
        "All_Zero": (
            all_zero_taxa.to_numpy()
        ),
    }
)


# ============================================================
# 16. EXTRACT DENOVO IDENTIFIER
# ============================================================

taxon_qc["Denovo_ID"] = (
    taxon_qc["Source_Taxon"]
    .str.extract(
        r"(d__denovo[^;]+)",
        expand=False
    )
)

missing_denovo = int(
    taxon_qc[
        "Denovo_ID"
    ].isna().sum()
)

duplicate_source_taxa = int(
    taxon_qc[
        "Source_Taxon"
    ].duplicated().sum()
)

duplicate_denovo_ids = int(
    taxon_qc[
        "Denovo_ID"
    ]
    .dropna()
    .duplicated()
    .sum()
)

print("\nTaxon identifier QC:")

print(
    f"Duplicate source taxon labels: "
    f"{duplicate_source_taxa}"
)

print(
    f"Missing d__denovo identifiers: "
    f"{missing_denovo}"
)

print(
    f"Duplicate d__denovo identifiers: "
    f"{duplicate_denovo_ids}"
)


# ============================================================
# 17. LIBRARY-SIZE OUTPUT TABLE
# ============================================================

library_size_df = pd.DataFrame(
    {
        "SampleID": (
            library_sizes.index
            .astype(str)
        ),
        "Library_Size": (
            library_sizes.values
        ),
    }
)

library_size_df = (
    library_size_df.merge(
        metadata[
            [
                "SampleID",
                "DiseaseState",
            ]
        ],
        on="SampleID",
        how="left",
        validate="one_to_one"
    )
)


# ============================================================
# 18. SAVE QC OUTPUTS
# ============================================================

threshold_file = (
    OUTPUT_DIR /
    "prevalence_threshold_comparison.csv"
)

taxon_qc_file = (
    OUTPUT_DIR /
    "raw_taxon_qc.csv"
)

library_file = (
    OUTPUT_DIR /
    "raw_library_sizes.csv"
)

metadata_qc_file = (
    OUTPUT_DIR /
    "raw_metadata_qc.csv"
)

summary_file = (
    OUTPUT_DIR /
    "raw_qc_summary.txt"
)


threshold_df.to_csv(
    threshold_file,
    index=False
)

taxon_qc.to_csv(
    taxon_qc_file,
    index=False
)

library_size_df.to_csv(
    library_file,
    index=False
)

metadata.to_csv(
    metadata_qc_file,
    index=False
)


# ============================================================
# 19. OVERALL PASS / FAIL
# ============================================================

critical_checks = {
    "No duplicate metadata SampleIDs":
        duplicate_metadata_ids == 0,

    "Expected groups only":
        len(unexpected_groups) == 0,

    "Both H and PAR are present":
        len(missing_expected_groups) == 0,

    "Sample sets match exactly":
        sample_sets_match,

    "No non-numeric counts":
        non_numeric_cells == 0,

    "No negative counts":
        negative_cells == 0,

    "No non-integer counts":
        non_integer_cells == 0,

    "No zero-library samples":
        len(zero_library_samples) == 0,
}


overall_pass = all(
    critical_checks.values()
)

print("\nCritical QC checks:")

for name, passed in critical_checks.items():

    print(
        f"{name}: "
        f"{'PASS' if passed else 'FAIL'}"
    )


# ============================================================
# 20. WRITE TEXT SUMMARY
# ============================================================

summary_lines = [
    "PAR_SCHEPERJANS RAW DATA QC",
    "=" * 64,
    "",
    f"Metadata samples: {metadata.shape[0]}",
    f"Raw taxa: {counts_numeric.shape[0]}",
    f"Count-table samples: {counts_numeric.shape[1]}",
    "",
    "GROUP COUNTS",
]

for group, count in group_counts.items():

    summary_lines.append(
        f"{group}: {count}"
    )


summary_lines.extend(
    [
        "",
        "SAMPLE ALIGNMENT",
        (
            f"Sample sets match exactly: "
            f"{sample_sets_match}"
        ),
        (
            f"Metadata-only samples: "
            f"{len(metadata_not_counts)}"
        ),
        (
            f"Count-only samples: "
            f"{len(counts_not_metadata)}"
        ),
        "",
        "COUNT QC",
        (
            f"Non-numeric cells: "
            f"{non_numeric_cells}"
        ),
        (
            f"Negative cells: "
            f"{negative_cells}"
        ),
        (
            f"Non-integer cells: "
            f"{non_integer_cells}"
        ),
        (
            f"Zero-library samples: "
            f"{len(zero_library_samples)}"
        ),
        "",
        "LIBRARY SIZE",
        (
            f"Mean: "
            f"{library_sizes.mean():.2f}"
        ),
        (
            f"Median: "
            f"{library_sizes.median():.2f}"
        ),
        (
            f"Minimum: "
            f"{library_sizes.min()}"
        ),
        (
            f"Maximum: "
            f"{library_sizes.max()}"
        ),
        "",
        "SPARSITY",
        (
            f"Zero fraction: "
            f"{zero_fraction:.6f}"
        ),
        (
            f"Zero percentage: "
            f"{zero_fraction * 100:.2f}%"
        ),
        "",
        "TAXON QC",
        (
            f"All-zero taxa: "
            f"{int(all_zero_taxa.sum())}"
        ),
        (
            f"Missing denovo identifiers: "
            f"{missing_denovo}"
        ),
        (
            f"Duplicate denovo identifiers: "
            f"{duplicate_denovo_ids}"
        ),
        "",
        "OVERALL RESULT",
        (
            "PASS"
            if overall_pass
            else "FAIL"
        ),
    ]
)

summary_file.write_text(
    "\n".join(summary_lines),
    encoding="utf-8"
)


# ============================================================
# 21. FINAL OUTPUT
# ============================================================

print("\n" + "=" * 72)

print(
    "OVERALL RAW-DATA AUDIT: "
    f"{'PASS' if overall_pass else 'FAIL'}"
)

print("=" * 72)

print("\nSaved QC files:")

print(threshold_file)
print(taxon_qc_file)
print(library_file)
print(metadata_qc_file)
print(summary_file)
