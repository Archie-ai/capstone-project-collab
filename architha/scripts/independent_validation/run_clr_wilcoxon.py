from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# ============================================================
# CONFIGURATION
# ============================================================

COUNTS_FILE = Path(
    "data/processed/simulation_null/null_counts.csv"
)

METADATA_FILE = Path(
    "data/processed/simulation_null/null_metadata.csv"
)

TRUTH_FILE = Path(
    "data/processed/simulation_null/null_ground_truth.csv"
)

OUTPUT_DIR = Path(
    "results/independent_validation/clr_wilcoxon"
)

PSEUDOCOUNT = 0.5
ALPHA = 0.05


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 1. LOAD INPUT DATA
# ============================================================

counts = pd.read_csv(
    COUNTS_FILE,
    index_col=0
)

metadata = pd.read_csv(
    METADATA_FILE
)

truth = pd.read_csv(
    TRUTH_FILE,
    index_col=0,
    keep_default_na=False
)


print("=" * 70)
print("CLR/WILCOXON INDEPENDENT VALIDATION")
print("=" * 70)

print("\nInput files:")
print(f"Counts:   {COUNTS_FILE}")
print(f"Metadata: {METADATA_FILE}")
print(f"Truth:    {TRUTH_FILE}")

print("\nInput dimensions:")
print(f"Counts:   {counts.shape}")
print(f"Metadata: {metadata.shape}")
print(f"Truth:    {truth.shape}")


# ============================================================
# 2. VALIDATE REQUIRED METADATA COLUMNS
# ============================================================

required_metadata_columns = {
    "SampleID",
    "Group"
}

missing_metadata_columns = (
    required_metadata_columns
    - set(metadata.columns)
)

if missing_metadata_columns:
    raise ValueError(
        "Missing required metadata columns: "
        f"{sorted(missing_metadata_columns)}"
    )


# ============================================================
# 3. VALIDATE SAMPLE IDENTIFIERS
# ============================================================

count_samples = list(
    counts.columns
)

metadata_samples = metadata[
    "SampleID"
].tolist()


if set(count_samples) != set(metadata_samples):

    missing_from_metadata = (
        set(count_samples)
        - set(metadata_samples)
    )

    missing_from_counts = (
        set(metadata_samples)
        - set(count_samples)
    )

    raise ValueError(
        "Sample IDs do not match.\n"
        f"Missing from metadata: "
        f"{sorted(missing_from_metadata)}\n"
        f"Missing from counts: "
        f"{sorted(missing_from_counts)}"
    )


# Reorder metadata so that its sample order exactly
# matches the columns of the count matrix.

metadata = (
    metadata
    .set_index("SampleID")
    .loc[count_samples]
    .reset_index()
)


samples_aligned = (
    list(counts.columns)
    == metadata["SampleID"].tolist()
)

print("\nSample alignment:")
print(
    "Counts and metadata aligned:",
    samples_aligned
)


if not samples_aligned:
    raise ValueError(
        "Counts and metadata are not aligned."
    )


# ============================================================
# 4. IDENTIFY CONTROL AND CASE SAMPLES
# ============================================================

control_samples = metadata.loc[
    metadata["Group"] == "Control",
    "SampleID"
].tolist()

case_samples = metadata.loc[
    metadata["Group"] == "Case",
    "SampleID"
].tolist()


print("\nGroup counts:")
print(f"Control: {len(control_samples)}")
print(f"Case:    {len(case_samples)}")


if len(control_samples) == 0:
    raise ValueError(
        "No Control samples detected."
    )

if len(case_samples) == 0:
    raise ValueError(
        "No Case samples detected."
    )


# ============================================================
# 5. VALIDATE RAW COUNTS
# ============================================================

all_numeric = all(
    pd.api.types.is_numeric_dtype(
        counts[column]
    )
    for column in counts.columns
)

if not all_numeric:
    raise ValueError(
        "Count matrix contains nonnumeric values."
    )


negative_count_total = int(
    (counts < 0)
    .sum()
    .sum()
)

if negative_count_total > 0:
    raise ValueError(
        f"Negative counts detected: "
        f"{negative_count_total}"
    )


integer_mask = np.isclose(
    counts.to_numpy(),
    np.round(
        counts.to_numpy()
    )
)

non_integer_total = int(
    (~integer_mask)
    .sum()
)

if non_integer_total > 0:
    raise ValueError(
        f"Non-integer counts detected: "
        f"{non_integer_total}"
    )


all_zero_taxa = (
    counts.sum(
        axis=1
    ) == 0
)


print("\nRaw-count QC:")
print(
    f"Total taxa: "
    f"{counts.shape[0]}"
)

print(
    "Negative counts:",
    negative_count_total
)

print(
    "Non-integer counts:",
    non_integer_total
)

print(
    "All-zero taxa:",
    int(
        all_zero_taxa.sum()
    )
)


# ============================================================
# 6. ADD PSEUDOCOUNT
# ============================================================

counts_with_pseudocount = (
    counts.astype(float)
    + PSEUDOCOUNT
)


print("\nCLR configuration:")
print(
    f"Pseudocount: "
    f"{PSEUDOCOUNT}"
)


# ============================================================
# 7. CLR TRANSFORMATION
# ============================================================

# Count table orientation:
#
# rows    = taxa
# columns = samples
#
# For each sample j:
#
# CLR_ij =
# log(x_ij + pseudocount)
# -
# mean_k[
#     log(x_kj + pseudocount)
# ]
#
# Therefore, the mean is calculated across taxa
# independently within each sample.

log_counts = np.log(
    counts_with_pseudocount
)


clr = log_counts.subtract(
    log_counts.mean(
        axis=0
    ),
    axis=1
)


# The CLR-transformed values for each sample should
# sum approximately to zero.

clr_column_sums = clr.sum(
    axis=0
)

maximum_abs_clr_sum = float(
    clr_column_sums
    .abs()
    .max()
)


print(
    "Maximum absolute CLR column sum:",
    f"{maximum_abs_clr_sum:.12g}"
)


if maximum_abs_clr_sum > 1e-8:
    print(
        "WARNING: CLR column sums are not "
        "numerically close to zero."
    )


clr.index.name = (
    counts.index.name
)


clr_file = (
    OUTPUT_DIR
    / "clr_transformed_matrix.csv"
)

clr.to_csv(
    clr_file
)


# ============================================================
# 8. RUN TAXON-BY-TAXON WILCOXON / MANN-WHITNEY TEST
# ============================================================

results_list = []


for taxon in clr.index:

    control_values = (
        clr.loc[
            taxon,
            control_samples
        ]
        .astype(float)
        .to_numpy()
    )

    case_values = (
        clr.loc[
            taxon,
            case_samples
        ]
        .astype(float)
        .to_numpy()
    )


    # Mann-Whitney U is the two-sample equivalent of
    # the Wilcoxon rank-sum test.
    #
    # The asymptotic implementation is used because
    # ties can occur in microbiome data.

    test_result = mannwhitneyu(
        case_values,
        control_values,
        alternative="two-sided",
        method="asymptotic"
    )


    control_mean_clr = float(
        np.mean(
            control_values
        )
    )

    case_mean_clr = float(
        np.mean(
            case_values
        )
    )

    control_median_clr = float(
        np.median(
            control_values
        )
    )

    case_median_clr = float(
        np.median(
            case_values
        )
    )


    effect_estimate = (
        case_mean_clr
        - control_mean_clr
    )


    if effect_estimate > 0:

        effect_direction = (
            "Case_higher"
        )

    elif effect_estimate < 0:

        effect_direction = (
            "Control_higher"
        )

    else:

        effect_direction = (
            "No_difference"
        )


    raw_taxon_counts = (
        counts.loc[
            taxon
        ]
    )


    present_samples = int(
        (
            raw_taxon_counts
            > 0
        ).sum()
    )


    present_control = int(
        (
            raw_taxon_counts[
                control_samples
            ]
            > 0
        ).sum()
    )


    present_case = int(
        (
            raw_taxon_counts[
                case_samples
            ]
            > 0
        ).sum()
    )


    results_list.append(
        {
            "Taxon": taxon,
            "Method": "CLR_Wilcoxon",
            "Tested": True,
            "P_Value": float(
                test_result.pvalue
            ),
            "Test_Statistic": float(
                test_result.statistic
            ),
            "Effect_Estimate": float(
                effect_estimate
            ),
            "Effect_Metric":
                "Case_minus_Control_mean_CLR",
            "Effect_Direction":
                effect_direction,
            "Control_Mean_CLR":
                control_mean_clr,
            "Case_Mean_CLR":
                case_mean_clr,
            "Control_Median_CLR":
                control_median_clr,
            "Case_Median_CLR":
                case_median_clr,
            "Present_Samples":
                present_samples,
            "Present_Control":
                present_control,
            "Present_Case":
                present_case,
            "All_Zero": bool(
                all_zero_taxa.loc[
                    taxon
                ]
            ),
            "Pseudocount":
                PSEUDOCOUNT
        }
    )


results = pd.DataFrame(
    results_list
)


# ============================================================
# 9. BENJAMINI-HOCHBERG MULTIPLE-TESTING CORRECTION
# ============================================================

reject, adjusted_p_values, _, _ = (
    multipletests(
        results[
            "P_Value"
        ],
        alpha=ALPHA,
        method="fdr_bh"
    )
)


results[
    "Adjusted_P_Value"
] = adjusted_p_values


results[
    "Significant"
] = reject


# ============================================================
# 10. PREPARE GROUND-TRUTH TABLE
# ============================================================

truth_reset = (
    truth
    .reset_index()
)


# Rename the first column Taxon if necessary.

first_truth_column = (
    truth_reset.columns[0]
)

if first_truth_column != "Taxon":

    truth_reset = (
        truth_reset.rename(
            columns={
                first_truth_column:
                    "Taxon"
            }
        )
    )


# Ensure taxon IDs are strings.

results[
    "Taxon"
] = (
    results[
        "Taxon"
    ]
    .astype(str)
)

truth_reset[
    "Taxon"
] = (
    truth_reset[
        "Taxon"
    ]
    .astype(str)
)


# ============================================================
# 11. MERGE RESULTS WITH KNOWN SIMULATION TRUTH
# ============================================================

results = results.merge(
    truth_reset,
    on="Taxon",
    how="left",
    validate="one_to_one"
)


# Confirm every analyzed taxon received truth information.

if results[
    "True_DA"
].isna().any():

    missing_truth_taxa = (
        results.loc[
            results[
                "True_DA"
            ].isna(),
            "Taxon"
        ]
        .tolist()
    )

    raise ValueError(
        "Ground truth missing for taxa: "
        f"{missing_truth_taxa}"
    )


# ============================================================
# 12. ADD ANALYSIS METADATA
# ============================================================

results[
    "Alpha"
] = ALPHA


results[
    "Scenario"
] = (
    "Null_Baseline_Low_Sparsity"
)


results[
    "Seed"
] = 42


# ============================================================
# 13. ORDER OUTPUT COLUMNS
# ============================================================

preferred_columns = [
    "Taxon",
    "Method",
    "Tested",
    "P_Value",
    "Adjusted_P_Value",
    "Significant",
    "Test_Statistic",
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
    "All_Zero",
    "Pseudocount",
    "Alpha",
    "True_DA",
    "True_Direction",
    "True_Effect_Size",
    "Scenario",
    "Seed"
]


existing_preferred_columns = [
    column
    for column in preferred_columns
    if column in results.columns
]


remaining_columns = [
    column
    for column in results.columns
    if column not in
    existing_preferred_columns
]


results = results[
    existing_preferred_columns
    + remaining_columns
]


# ============================================================
# 14. SORT RESULTS
# ============================================================

results = (
    results
    .sort_values(
        by=[
            "Adjusted_P_Value",
            "P_Value"
        ],
        ascending=[
            True,
            True
        ]
    )
    .reset_index(
        drop=True
    )
)


# ============================================================
# 15. SAVE FULL RESULTS
# ============================================================

results_file = (
    OUTPUT_DIR
    / "clr_wilcoxon_results.csv"
)


results.to_csv(
    results_file,
    index=False
)


# ============================================================
# 16. CALCULATE SUMMARY STATISTICS
# ============================================================

n_taxa_available = int(
    counts.shape[0]
)

n_tested = int(
    results[
        "Tested"
    ].sum()
)

n_significant = int(
    results[
        "Significant"
    ].sum()
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


significant_taxa = (
    results.loc[
        results[
            "Significant"
        ],
        "Taxon"
    ]
    .tolist()
)


# ============================================================
# 17. CALCULATE NULL PERFORMANCE
# ============================================================

true_da_count = int(
    pd.to_numeric(
        results[
            "True_DA"
        ],
        errors="coerce"
    )
    .fillna(0)
    .sum()
)


# Under the global null, every significant result
# is necessarily a false positive.

false_positive_count = (
    n_significant
)


# ============================================================
# 18. WRITE TEXT SUMMARY
# ============================================================

summary_lines = [
    "CLR/WILCOXON INDEPENDENT VALIDATION",
    "=" * 60,
    "",
    "INPUT",
    f"Count table: {COUNTS_FILE}",
    f"Metadata: {METADATA_FILE}",
    f"Truth: {TRUTH_FILE}",
    "",
    "CONFIGURATION",
    f"Pseudocount: {PSEUDOCOUNT}",
    "Transformation: CLR within each sample",
    (
        "Test: two-sided Mann-Whitney U "
        "(Wilcoxon rank-sum equivalent)"
    ),
    "SciPy method: asymptotic",
    (
        "Multiple-testing correction: "
        "Benjamini-Hochberg"
    ),
    f"Alpha: {ALPHA}",
    "",
    "DATASET",
    f"Taxa available: {n_taxa_available}",
    f"Control samples: {len(control_samples)}",
    f"Case samples: {len(case_samples)}",
    (
        "All-zero taxa: "
        f"{int(all_zero_taxa.sum())}"
    ),
    "",
    "RESULTS",
    f"Taxa tested: {n_tested}",
    f"Significant taxa: {n_significant}",
    (
        "Minimum raw p-value: "
        f"{minimum_raw_p:.10g}"
    ),
    (
        "Minimum adjusted p-value: "
        f"{minimum_adjusted_p:.10g}"
    ),
    "",
    "GROUND TRUTH",
    f"True DA taxa: {true_da_count}",
    (
        "False positives under null: "
        f"{false_positive_count}"
    ),
    "",
    "SIGNIFICANT TAXA",
    (
        ", ".join(
            significant_taxa
        )
        if significant_taxa
        else "None"
    ),
    "",
    "INTERPRETATION",
    (
        "Because this simulation contains no "
        "true differentially abundant taxa, "
        "any BH-significant taxon represents "
        "a false-positive discovery."
    )
]


summary_file = (
    OUTPUT_DIR
    / "clr_wilcoxon_summary.txt"
)


summary_file.write_text(
    "\n".join(
        summary_lines
    ),
    encoding="utf-8"
)


# ============================================================
# 19. PRINT RESULT SUMMARY
# ============================================================

print(
    "\n"
    + "=" * 70
)

print(
    "RESULT SUMMARY"
)

print(
    "=" * 70
)


print(
    f"Taxa available: "
    f"{n_taxa_available}"
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
    "Minimum raw p-value:",
    f"{minimum_raw_p:.10g}"
)

print(
    "Minimum adjusted p-value:",
    f"{minimum_adjusted_p:.10g}"
)

print(
    "True DA taxa:",
    true_da_count
)

print(
    "False positives under null:",
    false_positive_count
)


# ============================================================
# 20. PRINT TOP 15 TAXA
# ============================================================

print(
    "\nTop 15 taxa by adjusted p-value:"
)


display_columns = [
    "Taxon",
    "Present_Samples",
    "Present_Control",
    "Present_Case",
    "P_Value",
    "Adjusted_P_Value",
    "Effect_Estimate",
    "Effect_Direction",
    "Significant"
]


print(
    results[
        display_columns
    ]
    .head(15)
    .to_string(
        index=False
    )
)


# ============================================================
# 21. PRINT RAREST TAXA FOR COMPARISON
# ============================================================

print(
    "\nRarest taxa and their statistical results:"
)


rarest_taxa = (
    results
    .sort_values(
        by=[
            "Present_Samples",
            "Adjusted_P_Value"
        ]
    )
    .head(10)
)


print(
    rarest_taxa[
        display_columns
    ]
    .to_string(
        index=False
    )
)


# ============================================================
# 22. PRINT OUTPUT LOCATIONS
# ============================================================

print(
    "\nOutput files:"
)

print(
    clr_file
)

print(
    results_file
)

print(
    summary_file
)


print(
    "\nCLR/Wilcoxon independent "
    "validation completed."
)
