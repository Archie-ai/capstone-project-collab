#!/usr/bin/env python3

"""
CLR/Wilcoxon differential-abundance analysis
for the MicrobiomeHD par_scheperjans dataset.

Project:
Microbiome Abundance Consensus Machine

Canonical dataset:
691 prevalence-filtered taxa
148 samples
74 Control / 74 Parkinson's disease Case

Method:
1. Add pseudocount 0.5
2. Perform sample-wise centered log-ratio (CLR) transformation
3. Compare Case vs Control for every taxon using two-sided
   Mann-Whitney U / Wilcoxon rank-sum inference
4. Apply Benjamini-Hochberg correction
5. Significant threshold: adjusted p < 0.05

Effect:
Mean CLR(Case) - Mean CLR(Control)

Positive effect = Case_higher
Negative effect = Control_higher
"""

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# ============================================================
# 1. CONFIGURATION
# ============================================================

DATASET = "par_scheperjans"
PSEUDOCOUNT = 0.5
ALPHA = 0.05

COUNTS_FILE = Path(
    "data/processed/par_scheperjans/"
    "par_scheperjans_counts.csv"
)

METADATA_FILE = Path(
    "data/processed/par_scheperjans/"
    "par_scheperjans_metadata.csv"
)

TAXONOMY_FILE = Path(
    "data/processed/par_scheperjans/"
    "par_scheperjans_taxonomy_mapping.csv"
)

TAXON_QC_FILE = Path(
    "data/processed/par_scheperjans/"
    "par_scheperjans_taxon_qc.csv"
)

OUTPUT_DIR = Path(
    "results/real/par_scheperjans/"
    "clr_wilcoxon"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

RESULT_FILE = (
    OUTPUT_DIR /
    "clr_wilcoxon_results.csv"
)

SIGNIFICANT_FILE = (
    OUTPUT_DIR /
    "clr_wilcoxon_significant_taxa.csv"
)

CLR_MATRIX_FILE = (
    OUTPUT_DIR /
    "clr_transformed_matrix.csv"
)

SUMMARY_FILE = (
    OUTPUT_DIR /
    "clr_wilcoxon_summary.txt"
)


# ============================================================
# 2. LOAD CANONICAL DATA
# ============================================================

counts = pd.read_csv(
    COUNTS_FILE,
    index_col=0
)

metadata = pd.read_csv(
    METADATA_FILE
)

taxonomy = pd.read_csv(
    TAXONOMY_FILE
)

taxon_qc = pd.read_csv(
    TAXON_QC_FILE
)


print("=" * 72)
print("PAR_SCHEPERJANS CLR/WILCOXON ANALYSIS")
print("=" * 72)

print("\nCanonical input dimensions:")

print(
    f"Counts:   {counts.shape[0]} taxa x "
    f"{counts.shape[1]} samples"
)

print(
    f"Metadata: {metadata.shape[0]} samples x "
    f"{metadata.shape[1]} columns"
)

print(
    f"Taxonomy: {taxonomy.shape[0]} taxa x "
    f"{taxonomy.shape[1]} columns"
)

print(
    f"Taxon QC: {taxon_qc.shape[0]} taxa x "
    f"{taxon_qc.shape[1]} columns"
)


# ============================================================
# 3. VALIDATE INPUT STRUCTURE
# ============================================================

required_metadata_columns = {
    "SampleID",
    "Group",
}

missing_metadata_columns = (
    required_metadata_columns -
    set(metadata.columns)
)

if missing_metadata_columns:
    raise ValueError(
        "Missing metadata columns: "
        f"{sorted(missing_metadata_columns)}"
    )


metadata["SampleID"] = (
    metadata["SampleID"]
    .astype(str)
    .str.strip()
)

metadata["Group"] = (
    metadata["Group"]
    .astype(str)
    .str.strip()
)

counts.columns = (
    counts.columns
    .astype(str)
    .str.strip()
)

counts.index = (
    counts.index
    .astype(str)
    .str.strip()
)


if metadata["SampleID"].duplicated().any():
    raise ValueError(
        "Duplicate metadata SampleIDs detected."
    )

if counts.index.duplicated().any():
    raise ValueError(
        "Duplicate Taxon IDs detected."
    )

if counts.columns.duplicated().any():
    raise ValueError(
        "Duplicate count-table sample IDs detected."
    )


metadata_sample_set = set(
    metadata["SampleID"]
)

count_sample_set = set(
    counts.columns
)

if metadata_sample_set != count_sample_set:

    metadata_only = sorted(
        metadata_sample_set -
        count_sample_set
    )

    counts_only = sorted(
        count_sample_set -
        metadata_sample_set
    )

    raise ValueError(
        "Count-table and metadata sample sets do not match.\n"
        f"Metadata-only: {metadata_only}\n"
        f"Counts-only: {counts_only}"
    )


# ============================================================
# 4. ALIGN METADATA TO COUNT COLUMNS
#
# IMPORTANT:
# keep SampleID as an explicit column using drop=False.
# This avoids the reset_index naming issue encountered earlier.
# ============================================================

metadata = (
    metadata
    .set_index(
        "SampleID",
        drop=False
    )
    .reindex(
        counts.columns
    )
    .copy()
)


if metadata["SampleID"].isna().any():
    raise ValueError(
        "Metadata reindexing produced missing SampleIDs."
    )


sample_alignment = (
    list(counts.columns)
    ==
    metadata["SampleID"].tolist()
)


print("\nSample alignment:")

print(
    f"Counts and metadata aligned: "
    f"{sample_alignment}"
)

print(
    f"Metadata SampleID column retained: "
    f"{'SampleID' in metadata.columns}"
)


if not sample_alignment:
    raise ValueError(
        "Sample alignment failed."
    )


# ============================================================
# 5. VALIDATE GROUPS
# ============================================================

observed_groups = set(
    metadata["Group"]
)

expected_groups = {
    "Control",
    "Case",
}


if observed_groups != expected_groups:
    raise ValueError(
        "Unexpected group labels.\n"
        f"Observed: {sorted(observed_groups)}\n"
        f"Expected: {sorted(expected_groups)}"
    )


group_counts = (
    metadata["Group"]
    .value_counts()
)


n_control = int(
    group_counts.get(
        "Control",
        0
    )
)

n_case = int(
    group_counts.get(
        "Case",
        0
    )
)


print("\nGroup counts:")

print(
    f"Control: {n_control}"
)

print(
    f"Case:    {n_case}"
)


if n_control != 74:
    raise ValueError(
        f"Expected 74 Control samples; found {n_control}."
    )

if n_case != 74:
    raise ValueError(
        f"Expected 74 Case samples; found {n_case}."
    )


# ============================================================
# 6. VALIDATE COUNTS
# ============================================================

counts_numeric = counts.apply(
    pd.to_numeric,
    errors="coerce"
)


missing_or_non_numeric = int(
    counts_numeric
    .isna()
    .sum()
    .sum()
)


values = counts_numeric.to_numpy(
    dtype=float
)


non_finite_count = int(
    (~np.isfinite(values)).sum()
)


negative_count = int(
    (
        counts_numeric < 0
    )
    .sum()
    .sum()
)


non_integer_count = int(
    (
        ~np.isclose(
            values,
            np.round(values)
        )
    ).sum()
)


all_zero_taxa = int(
    (
        counts_numeric.sum(
            axis=1
        ) == 0
    ).sum()
)


zero_library_samples = int(
    (
        counts_numeric.sum(
            axis=0
        ) == 0
    ).sum()
)


print("\nCount QC:")

print(
    f"Missing/non-numeric cells: "
    f"{missing_or_non_numeric}"
)

print(
    f"Non-finite cells: "
    f"{non_finite_count}"
)

print(
    f"Negative cells: "
    f"{negative_count}"
)

print(
    f"Non-integer cells: "
    f"{non_integer_count}"
)

print(
    f"All-zero taxa: "
    f"{all_zero_taxa}"
)

print(
    f"Zero-library samples: "
    f"{zero_library_samples}"
)


if missing_or_non_numeric != 0:
    raise ValueError(
        "Missing or non-numeric counts detected."
    )

if non_finite_count != 0:
    raise ValueError(
        "Non-finite counts detected."
    )

if negative_count != 0:
    raise ValueError(
        "Negative counts detected."
    )

if non_integer_count != 0:
    raise ValueError(
        "Non-integer counts detected."
    )

if all_zero_taxa != 0:
    raise ValueError(
        "All-zero taxa detected."
    )

if zero_library_samples != 0:
    raise ValueError(
        "Zero-library samples detected."
    )


counts_numeric = (
    counts_numeric
    .round()
    .astype("int64")
)


# ============================================================
# 7. IDENTIFY GROUP SAMPLE IDS
# ============================================================

control_ids = (
    metadata.loc[
        metadata["Group"]
        ==
        "Control",
        "SampleID"
    ]
    .tolist()
)

case_ids = (
    metadata.loc[
        metadata["Group"]
        ==
        "Case",
        "SampleID"
    ]
    .tolist()
)


if len(control_ids) != 74:
    raise ValueError(
        "Control sample-ID extraction failed."
    )

if len(case_ids) != 74:
    raise ValueError(
        "Case sample-ID extraction failed."
    )


# ============================================================
# 8. CLR TRANSFORMATION
#
# CLR(x_i) = ln(x_i) - mean(ln(x))
#
# Calculated independently within each sample.
# ============================================================

counts_plus = (
    counts_numeric.astype(float)
    +
    PSEUDOCOUNT
)


log_counts = np.log(
    counts_plus
)


clr_matrix = (
    log_counts
    -
    log_counts.mean(
        axis=0
    )
)


clr_matrix.index.name = "Taxon"


# ============================================================
# 9. CLR SANITY CHECK
# ============================================================

clr_column_sums = (
    clr_matrix.sum(
        axis=0
    )
)


maximum_abs_clr_sum = float(
    clr_column_sums
    .abs()
    .max()
)


print("\nCLR transformation:")

print(
    f"Pseudocount: "
    f"{PSEUDOCOUNT}"
)

print(
    "Maximum absolute CLR column sum: "
    f"{maximum_abs_clr_sum:.12g}"
)


if maximum_abs_clr_sum > 1e-8:
    raise RuntimeError(
        "CLR column-sum sanity check failed."
    )


# ============================================================
# 10. PREVALENCE INFORMATION
# ============================================================

present_total = (
    counts_numeric > 0
).sum(
    axis=1
)

present_control = (
    counts_numeric[
        control_ids
    ] > 0
).sum(
    axis=1
)

present_case = (
    counts_numeric[
        case_ids
    ] > 0
).sum(
    axis=1
)


# ============================================================
# 11. WILCOXON / MANN-WHITNEY TEST
# ============================================================

rows = []


for taxon in clr_matrix.index:

    control_values = (
        clr_matrix.loc[
            taxon,
            control_ids
        ]
        .to_numpy(
            dtype=float
        )
    )

    case_values = (
        clr_matrix.loc[
            taxon,
            case_ids
        ]
        .to_numpy(
            dtype=float
        )
    )


    test = mannwhitneyu(
        case_values,
        control_values,
        alternative="two-sided",
        method="asymptotic"
    )


    control_mean = float(
        np.mean(
            control_values
        )
    )

    case_mean = float(
        np.mean(
            case_values
        )
    )


    control_median = float(
        np.median(
            control_values
        )
    )

    case_median = float(
        np.median(
            case_values
        )
    )


    effect = (
        case_mean -
        control_mean
    )


    if effect > 0:

        direction = (
            "Case_higher"
        )

    elif effect < 0:

        direction = (
            "Control_higher"
        )

    else:

        direction = (
            "No_difference"
        )


    rows.append(
        {
            "Taxon":
                taxon,

            "Method":
                "CLR_Wilcoxon",

            "Tested":
                True,

            "U_Statistic":
                float(
                    test.statistic
                ),

            "P_Value":
                float(
                    test.pvalue
                ),

            "Control_Mean_CLR":
                control_mean,

            "Case_Mean_CLR":
                case_mean,

            "Control_Median_CLR":
                control_median,

            "Case_Median_CLR":
                case_median,

            "Effect_Estimate":
                effect,

            "Effect_Metric":
                (
                    "CLR mean difference "
                    "(Case-Control)"
                ),

            "Effect_Direction":
                direction,

            "Present_Samples":
                int(
                    present_total.loc[
                        taxon
                    ]
                ),

            "Present_Control":
                int(
                    present_control.loc[
                        taxon
                    ]
                ),

            "Present_Case":
                int(
                    present_case.loc[
                        taxon
                    ]
                ),
        }
    )


results = pd.DataFrame(
    rows
)


# ============================================================
# 12. BENJAMINI-HOCHBERG CORRECTION
# ============================================================

reject, adjusted_p, _, _ = (
    multipletests(
        results["P_Value"],
        alpha=ALPHA,
        method="fdr_bh"
    )
)


results[
    "Adjusted_P_Value"
] = adjusted_p


results[
    "Significant"
] = reject


results[
    "Alpha"
] = ALPHA


results[
    "P_Adjustment"
] = (
    "Benjamini-Hochberg"
)


results[
    "Dataset"
] = DATASET


results[
    "Scenario"
] = "Real_Data"


# ============================================================
# 13. ADD TAXONOMY
# ============================================================

required_taxonomy_columns = {
    "Taxon",
}

missing_taxonomy_columns = (
    required_taxonomy_columns -
    set(taxonomy.columns)
)

if missing_taxonomy_columns:
    raise ValueError(
        "Taxonomy mapping is missing Taxon."
    )


taxonomy_columns = [
    column
    for column in [
        "Taxon",
        "Denovo_ID",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Source_Taxon",
    ]
    if column in taxonomy.columns
]


results = results.merge(
    taxonomy[
        taxonomy_columns
    ],
    on="Taxon",
    how="left",
    validate="one_to_one"
)


missing_taxonomy_matches = int(
    results[
        "Denovo_ID"
    ].isna().sum()
) if "Denovo_ID" in results.columns else 0


print("\nTaxonomy merge:")

print(
    f"Missing taxonomy matches: "
    f"{missing_taxonomy_matches}"
)


if missing_taxonomy_matches != 0:
    raise ValueError(
        "Some analysis taxa failed taxonomy mapping."
    )


# ============================================================
# 14. FINAL RESULT COLUMN ORDER
# ============================================================

preferred_columns = [
    "Taxon",
    "Method",
    "Tested",
    "U_Statistic",
    "P_Value",
    "Adjusted_P_Value",
    "Significant",
    "Effect_Estimate",
    "Effect_Metric",
    "Effect_Direction",
    "Control_Mean_CLR",
    "Case_Mean_CLR",
    "Control_Median_CLR",
    "Case_Median_CLR",
    "Present_Samples",
    "Present_Control",
    "Present_Case",
    "Alpha",
    "P_Adjustment",
    "Dataset",
    "Scenario",
    "Denovo_ID",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Source_Taxon",
]


results = results[
    [
        column
        for column in preferred_columns
        if column in results.columns
    ]
]


# ============================================================
# 15. SORT RESULTS
# ============================================================

results = (
    results
    .sort_values(
        by=[
            "Adjusted_P_Value",
            "P_Value",
            "Taxon",
        ],
        ascending=[
            True,
            True,
            True,
        ]
    )
    .reset_index(
        drop=True
    )
)


significant_results = (
    results.loc[
        results[
            "Significant"
        ]
    ]
    .copy()
)


# ============================================================
# 16. RESULT SUMMARY
# ============================================================

n_available = int(
    len(results)
)

n_tested = int(
    results["Tested"].sum()
)

n_significant = int(
    results["Significant"].sum()
)


minimum_raw_p = float(
    results[
        "P_Value"
    ].min()
)


minimum_adjusted_p = float(
    results[
        "Adjusted_P_Value"
    ].min()
)


case_higher_sig = int(
    (
        significant_results[
            "Effect_Direction"
        ]
        ==
        "Case_higher"
    ).sum()
)


control_higher_sig = int(
    (
        significant_results[
            "Effect_Direction"
        ]
        ==
        "Control_higher"
    ).sum()
)


print("\n" + "=" * 72)
print("RESULT SUMMARY")
print("=" * 72)


print(
    f"Taxa available: "
    f"{n_available}"
)

print(
    f"Taxa tested: "
    f"{n_tested}"
)

print(
    f"Significant taxa at BH q < "
    f"{ALPHA}: "
    f"{n_significant}"
)

print(
    f"Case-higher significant taxa: "
    f"{case_higher_sig}"
)

print(
    f"Control-higher significant taxa: "
    f"{control_higher_sig}"
)

print(
    f"Minimum raw p-value: "
    f"{minimum_raw_p:.10g}"
)

print(
    f"Minimum adjusted p-value: "
    f"{minimum_adjusted_p:.10g}"
)


# ============================================================
# 17. DISPLAY TOP 20
# ============================================================

print(
    "\nTop 20 taxa by adjusted p-value:"
)


top_columns = [
    "Taxon",
    "Denovo_ID",
    "Present_Samples",
    "Present_Control",
    "Present_Case",
    "Effect_Estimate",
    "Effect_Direction",
    "P_Value",
    "Adjusted_P_Value",
    "Significant",
    "Phylum",
    "Family",
    "Genus",
]


top_columns = [
    column
    for column in top_columns
    if column in results.columns
]


print(
    results[
        top_columns
    ]
    .head(20)
    .to_string(
        index=False
    )
)


# ============================================================
# 18. SAVE OUTPUTS
# ============================================================

results.to_csv(
    RESULT_FILE,
    index=False
)

significant_results.to_csv(
    SIGNIFICANT_FILE,
    index=False
)

clr_matrix.to_csv(
    CLR_MATRIX_FILE
)


# ============================================================
# 19. WRITE SUMMARY
# ============================================================

summary_lines = [
    "PAR_SCHEPERJANS CLR/WILCOXON SUMMARY",
    "=" * 68,
    "",
    "DATASET",
    (
        f"Taxa analyzed: "
        f"{n_available}"
    ),
    (
        f"Control samples: "
        f"{n_control}"
    ),
    (
        f"Case samples: "
        f"{n_case}"
    ),
    "",
    "METHOD",
    (
        f"CLR pseudocount: "
        f"{PSEUDOCOUNT}"
    ),
    (
        "Test: two-sided "
        "Mann-Whitney U / "
        "Wilcoxon rank-sum"
    ),
    (
        "SciPy implementation: "
        "asymptotic method"
    ),
    (
        "Effect: mean CLR(Case) - "
        "mean CLR(Control)"
    ),
    (
        "Multiple-testing correction: "
        "Benjamini-Hochberg"
    ),
    (
        f"Alpha: "
        f"{ALPHA}"
    ),
    "",
    "CLR QC",
    (
        "Maximum absolute sample-wise "
        "CLR sum: "
        f"{maximum_abs_clr_sum:.12g}"
    ),
    "",
    "RESULTS",
    (
        f"Taxa tested: "
        f"{n_tested}"
    ),
    (
        f"Significant taxa: "
        f"{n_significant}"
    ),
    (
        f"Case-higher significant taxa: "
        f"{case_higher_sig}"
    ),
    (
        "Control-higher significant taxa: "
        f"{control_higher_sig}"
    ),
    (
        f"Minimum raw p-value: "
        f"{minimum_raw_p:.10g}"
    ),
    (
        f"Minimum adjusted p-value: "
        f"{minimum_adjusted_p:.10g}"
    ),
    "",
    "INTERPRETATION NOTE",
    (
        "This is a real microbiome dataset "
        "with unknown ground truth."
    ),
    (
        "Significant taxa represent statistical "
        "discoveries rather than confirmed "
        "biological truth."
    ),
    (
        "Results will be compared with ALDEx2, "
        "ANCOM-BC2, and DESeq2 for concordance "
        "and consensus analysis."
    ),
]


SUMMARY_FILE.write_text(
    "\n".join(
        summary_lines
    ),
    encoding="utf-8"
)


# ============================================================
# 20. COMPLETE
# ============================================================

print("\nOutput files:")

print(
    RESULT_FILE
)

print(
    SIGNIFICANT_FILE
)

print(
    CLR_MATRIX_FILE
)

print(
    SUMMARY_FILE
)


print("\n" + "=" * 72)

print(
    "CLR/WILCOXON REAL-DATA ANALYSIS: PASS"
)

print("=" * 72)
