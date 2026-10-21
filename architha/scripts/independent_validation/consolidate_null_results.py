#!/usr/bin/env python3

"""
Consolidate independent null-validation results from:

1. CLR/Wilcoxon
2. ALDEx2
3. ANCOM-BC2
4. DESeq2

Project:
Microbiome Abundance Consensus Machine

Purpose:
Create standardized long-format and wide-format result tables from
the same validated null dataset.

The simulation contains zero true differentially abundant taxa.
Therefore, every significant discovery in this dataset is a known
false-positive discovery.

Primary method conventions used here:
- CLR/Wilcoxon: BH-adjusted Wilcoxon result
- ALDEx2: Wilcoxon expected BH result (wi.eBH equivalent)
  Welch results are retained as diagnostics
- ANCOM-BC2: robust significance call when available
- DESeq2: primary result with standard independent-filtering behavior

The script preserves method-specific effect metrics rather than
assuming that all effect estimates are directly comparable.
"""

from pathlib import Path
import itertools
import numpy as np
import pandas as pd


# ============================================================
# 1. CONFIGURATION
# ============================================================

ALPHA = 0.05

CLR_FILE = Path(
    "results/independent_validation/"
    "clr_wilcoxon/clr_wilcoxon_results.csv"
)

ALDEX_FILE = Path(
    "results/independent_validation/"
    "aldex2/aldex2_results.csv"
)

ANCOM_FILE = Path(
    "results/independent_validation/"
    "ancombc2/ancombc2_standardized_results.csv"
)

DESEQ_FILE = Path(
    "results/independent_validation/"
    "deseq2/deseq2_standardized_results.csv"
)

TRUTH_FILE = Path(
    "data/processed/simulation_null/"
    "null_ground_truth.csv"
)

OUTPUT_DIR = Path(
    "results/independent_validation/consolidated_null"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. OUTPUT FILES
# ============================================================

LONG_FILE = OUTPUT_DIR / "null_four_method_long.csv"

WIDE_FILE = OUTPUT_DIR / "null_four_method_comparison.csv"

METHOD_SUMMARY_FILE = (
    OUTPUT_DIR / "null_method_summary.csv"
)

CONSENSUS_FILE = (
    OUTPUT_DIR / "null_consensus_all_taxa.csv"
)

TAXON033_FILE = (
    OUTPUT_DIR / "null_taxon_033_comparison.csv"
)

JACCARD_FILE = (
    OUTPUT_DIR / "null_significant_jaccard_matrix.csv"
)

CALL_AGREEMENT_FILE = (
    OUTPUT_DIR / "null_significance_call_agreement_matrix.csv"
)

SUMMARY_FILE = (
    OUTPUT_DIR / "null_consolidation_summary.txt"
)


# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

def require_file(path: Path) -> None:
    """Stop immediately if a required input file is missing."""

    if not path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{path}"
        )


def find_column(
    df: pd.DataFrame,
    candidates,
    required=True
):
    """
    Return the first matching column.

    Matching is attempted:
    1. exactly;
    2. case-insensitively.
    """

    for candidate in candidates:
        if candidate in df.columns:
            return candidate

    lower_lookup = {
        str(column).lower(): column
        for column in df.columns
    }

    for candidate in candidates:
        key = str(candidate).lower()

        if key in lower_lookup:
            return lower_lookup[key]

    if required:
        raise KeyError(
            "Could not find any of these columns:\n"
            f"{candidates}\n\n"
            "Available columns:\n"
            f"{list(df.columns)}"
        )

    return None


def convert_bool(series: pd.Series) -> pd.Series:
    """
    Convert common boolean representations to nullable Boolean dtype.
    """

    if pd.api.types.is_bool_dtype(series):
        return series.astype("boolean")

    if pd.api.types.is_numeric_dtype(series):
        return series.map(
            lambda x: pd.NA if pd.isna(x) else bool(x)
        ).astype("boolean")

    mapping = {
        "true": True,
        "false": False,
        "t": True,
        "f": False,
        "1": True,
        "0": False,
        "yes": True,
        "no": False,
    }

    return (
        series.astype("string")
        .str.strip()
        .str.lower()
        .map(mapping)
        .astype("boolean")
    )


def numeric_series(
    df: pd.DataFrame,
    column_name
) -> pd.Series:
    """
    Safely return numeric values or an all-NA numeric series.
    """

    if column_name is None:
        return pd.Series(
            np.nan,
            index=df.index,
            dtype=float
        )

    return pd.to_numeric(
        df[column_name],
        errors="coerce"
    )


def boolean_series(
    df: pd.DataFrame,
    column_name
) -> pd.Series:
    """
    Safely return nullable boolean values or all NA.
    """

    if column_name is None:
        return pd.Series(
            pd.NA,
            index=df.index,
            dtype="boolean"
        )

    return convert_bool(
        df[column_name]
    )


def derive_significance(
    q_values: pd.Series,
    supplied_significance: pd.Series | None = None
) -> pd.Series:
    """
    Prefer a supplied significance flag if available.
    Otherwise calculate q < ALPHA.
    """

    if supplied_significance is not None:

        result = supplied_significance.copy()

        missing = result.isna()

        result.loc[missing] = (
            q_values.loc[missing] < ALPHA
        )

        return result.astype("boolean")

    return (
        q_values < ALPHA
    ).astype("boolean")


def get_taxon_column(
    df: pd.DataFrame
) -> str:
    """Identify the taxon/feature identifier column."""

    candidates = [
        "Taxon",
        "taxon",
        "Feature",
        "feature",
        "FeatureID",
        "OTU",
    ]

    column = find_column(
        df,
        candidates,
        required=False
    )

    if column is not None:
        return column

    first_column = df.columns[0]

    if str(first_column).startswith("Unnamed"):
        return first_column

    raise KeyError(
        "Could not identify the taxon column.\n"
        f"Available columns: {list(df.columns)}"
    )


def normalize_taxon_column(
    df: pd.DataFrame
) -> pd.DataFrame:
    """Rename whichever identifier column exists to Taxon."""

    df = df.copy()

    taxon_column = get_taxon_column(df)

    if taxon_column != "Taxon":
        df = df.rename(
            columns={
                taxon_column: "Taxon"
            }
        )

    df["Taxon"] = (
        df["Taxon"]
        .astype(str)
        .str.strip()
    )

    return df


def safe_min(series: pd.Series):
    """Minimum finite numeric value; otherwise NaN."""

    values = pd.to_numeric(
        series,
        errors="coerce"
    )

    values = values[
        np.isfinite(values)
    ]

    if len(values) == 0:
        return np.nan

    return float(
        values.min()
    )


# ============================================================
# 4. VERIFY REQUIRED FILES
# ============================================================

required_files = [
    CLR_FILE,
    ALDEX_FILE,
    ANCOM_FILE,
    DESEQ_FILE,
    TRUTH_FILE,
]

for required_file in required_files:
    require_file(required_file)


# ============================================================
# 5. LOAD RAW RESULT TABLES
# ============================================================

clr_raw = normalize_taxon_column(
    pd.read_csv(CLR_FILE)
)

aldex_raw = normalize_taxon_column(
    pd.read_csv(ALDEX_FILE)
)

ancom_raw = normalize_taxon_column(
    pd.read_csv(ANCOM_FILE)
)

deseq_raw = normalize_taxon_column(
    pd.read_csv(DESEQ_FILE)
)

truth_raw = normalize_taxon_column(
    pd.read_csv(
        TRUTH_FILE,
        keep_default_na=False
    )
)


print("=" * 72)
print("CONSOLIDATED NULL COMPARISON")
print("=" * 72)

print("\nInput result dimensions:")

print(
    f"CLR/Wilcoxon: {clr_raw.shape}"
)

print(
    f"ALDEx2:       {aldex_raw.shape}"
)

print(
    f"ANCOM-BC2:    {ancom_raw.shape}"
)

print(
    f"DESeq2:       {deseq_raw.shape}"
)

print(
    f"Truth:        {truth_raw.shape}"
)


# ============================================================
# 6. VALIDATE TAXON IDENTIFIERS
# ============================================================

truth_taxa = set(
    truth_raw["Taxon"]
)

method_taxa = {
    "CLR_Wilcoxon": set(
        clr_raw["Taxon"]
    ),
    "ALDEx2": set(
        aldex_raw["Taxon"]
    ),
    "ANCOM_BC2": set(
        ancom_raw["Taxon"]
    ),
    "DESeq2": set(
        deseq_raw["Taxon"]
    ),
}

print("\nTaxon-universe validation:")

for method, taxa in method_taxa.items():

    missing = truth_taxa - taxa
    extra = taxa - truth_taxa

    print(
        f"{method}: "
        f"{len(taxa)} taxa | "
        f"missing from method={len(missing)} | "
        f"extra={len(extra)}"
    )

    if missing:
        raise ValueError(
            f"{method} is missing truth taxa:\n"
            f"{sorted(missing)}"
        )

    if extra:
        raise ValueError(
            f"{method} contains unexpected taxa:\n"
            f"{sorted(extra)}"
        )


if len(truth_taxa) != 100:

    raise ValueError(
        "Expected exactly 100 taxa in the current null truth table; "
        f"found {len(truth_taxa)}."
    )


# ============================================================
# 7. STANDARDIZE CLR/WILCOXON
# ============================================================

clr_p_col = find_column(
    clr_raw,
    [
        "P_Value",
        "p_value",
        "PValue",
        "pvalue",
        "Raw_P_Value",
    ]
)

clr_q_col = find_column(
    clr_raw,
    [
        "Adjusted_P_Value",
        "Adjusted_P",
        "BH_Adjusted_P",
        "q_value",
        "padj",
    ]
)

clr_sig_col = find_column(
    clr_raw,
    [
        "Significant",
        "BH_Significant",
        "Primary_Significant",
    ],
    required=False
)

clr_tested_col = find_column(
    clr_raw,
    [
        "Tested",
    ],
    required=False
)

clr_effect_col = find_column(
    clr_raw,
    [
        "Effect_Estimate",
        "Case_minus_Control_CLR",
        "CLR_Effect",
        "Effect",
    ],
    required=False
)

clr_direction_col = find_column(
    clr_raw,
    [
        "Effect_Direction",
        "Direction",
    ],
    required=False
)


clr_p = numeric_series(
    clr_raw,
    clr_p_col
)

clr_q = numeric_series(
    clr_raw,
    clr_q_col
)

if clr_sig_col is not None:

    clr_supplied_sig = boolean_series(
        clr_raw,
        clr_sig_col
    )

else:

    clr_supplied_sig = None


clr_sig = derive_significance(
    clr_q,
    clr_supplied_sig
)


if clr_tested_col is not None:

    clr_tested = boolean_series(
        clr_raw,
        clr_tested_col
    )

else:

    clr_tested = (
        clr_p.notna()
    ).astype("boolean")


clr_effect = numeric_series(
    clr_raw,
    clr_effect_col
)


if clr_direction_col is not None:

    clr_direction = (
        clr_raw[clr_direction_col]
        .astype("string")
    )

else:

    clr_direction = pd.Series(
        pd.NA,
        index=clr_raw.index,
        dtype="string"
    )


clr_long = pd.DataFrame(
    {
        "Taxon": clr_raw["Taxon"],
        "Method": "CLR_Wilcoxon",
        "Primary_Statistic": "Wilcoxon_rank_sum",
        "Tested": clr_tested,
        "P_Value": clr_p,
        "Adjusted_P_Value": clr_q,
        "Significant": clr_sig,
        "Effect_Estimate": clr_effect,
        "Effect_Metric": (
            "CLR mean difference (Case-Control)"
        ),
        "Effect_Direction": clr_direction,
    }
)


# ============================================================
# 8. STANDARDIZE ALDEX2
#
# Primary field for the current four-method comparison:
# Wilcoxon expected p/BH values.
#
# Welch values are retained later in the wide table.
# ============================================================

aldex_p_col = find_column(
    aldex_raw,
    [
        "Wilcoxon_P_Value",
        "wi.ep",
    ]
)

aldex_q_col = find_column(
    aldex_raw,
    [
        "Wilcoxon_Adjusted_P_Value",
        "wi.eBH",
    ]
)

aldex_sig_col = find_column(
    aldex_raw,
    [
        "Wilcoxon_Significant",
    ],
    required=False
)

aldex_effect_col = find_column(
    aldex_raw,
    [
        "Case_minus_Control_CLR",
        "ALDEx2_Effect",
        "Effect_Estimate",
    ],
    required=False
)

aldex_direction_col = find_column(
    aldex_raw,
    [
        "Effect_Direction",
        "Direction",
    ],
    required=False
)

aldex_p = numeric_series(
    aldex_raw,
    aldex_p_col
)

aldex_q = numeric_series(
    aldex_raw,
    aldex_q_col
)

if aldex_sig_col is not None:

    aldex_supplied_sig = boolean_series(
        aldex_raw,
        aldex_sig_col
    )

else:

    aldex_supplied_sig = None


aldex_sig = derive_significance(
    aldex_q,
    aldex_supplied_sig
)

aldex_tested = (
    aldex_p.notna()
).astype("boolean")

aldex_effect = numeric_series(
    aldex_raw,
    aldex_effect_col
)


if aldex_direction_col is not None:

    aldex_direction = (
        aldex_raw[aldex_direction_col]
        .astype("string")
    )

else:

    aldex_direction = pd.Series(
        pd.NA,
        index=aldex_raw.index,
        dtype="string"
    )


if (
    aldex_effect_col ==
    "Case_minus_Control_CLR"
):

    aldex_effect_metric = (
        "ALDEx2 CLR difference (Case-Control)"
    )

else:

    aldex_effect_metric = (
        "ALDEx2 native effect statistic"
    )


aldex_long = pd.DataFrame(
    {
        "Taxon": aldex_raw["Taxon"],
        "Method": "ALDEx2",
        "Primary_Statistic": "Wilcoxon_expected",
        "Tested": aldex_tested,
        "P_Value": aldex_p,
        "Adjusted_P_Value": aldex_q,
        "Significant": aldex_sig,
        "Effect_Estimate": aldex_effect,
        "Effect_Metric": aldex_effect_metric,
        "Effect_Direction": aldex_direction,
    }
)


# ============================================================
# 9. STANDARDIZE ANCOM-BC2
#
# Robust_Significant is used as the primary final call when
# available because pseudo-count sensitivity analysis was run.
#
# Standard BH significance is retained in the wide table.
# ============================================================

ancom_p_col = find_column(
    ancom_raw,
    [
        "P_Value",
        "p_GroupCase",
    ]
)

ancom_q_col = find_column(
    ancom_raw,
    [
        "Adjusted_P_Value",
        "q_GroupCase",
    ]
)

ancom_standard_sig_col = find_column(
    ancom_raw,
    [
        "Standard_Significant",
        "Package_Diff_Call",
    ],
    required=False
)

ancom_robust_sig_col = find_column(
    ancom_raw,
    [
        "Robust_Significant",
        "diff_robust_GroupCase",
    ],
    required=False
)

ancom_effect_col = find_column(
    ancom_raw,
    [
        "Log_Fold_Change",
        "lfc_GroupCase",
    ],
    required=False
)

ancom_direction_col = find_column(
    ancom_raw,
    [
        "Effect_Direction",
        "Direction",
    ],
    required=False
)

ancom_tested_col = find_column(
    ancom_raw,
    [
        "Tested",
    ],
    required=False
)


ancom_p = numeric_series(
    ancom_raw,
    ancom_p_col
)

ancom_q = numeric_series(
    ancom_raw,
    ancom_q_col
)


if ancom_standard_sig_col is not None:

    ancom_standard_sig = boolean_series(
        ancom_raw,
        ancom_standard_sig_col
    )

else:

    ancom_standard_sig = (
        ancom_q < ALPHA
    ).astype("boolean")


if ancom_robust_sig_col is not None:

    ancom_primary_sig = boolean_series(
        ancom_raw,
        ancom_robust_sig_col
    )

else:

    ancom_primary_sig = ancom_standard_sig.copy()


if ancom_tested_col is not None:

    ancom_tested = boolean_series(
        ancom_raw,
        ancom_tested_col
    )

else:

    ancom_tested = (
        ancom_p.notna()
    ).astype("boolean")


ancom_effect = numeric_series(
    ancom_raw,
    ancom_effect_col
)


if ancom_direction_col is not None:

    ancom_direction = (
        ancom_raw[ancom_direction_col]
        .astype("string")
    )

else:

    ancom_direction = pd.Series(
        pd.NA,
        index=ancom_raw.index,
        dtype="string"
    )


ancom_long = pd.DataFrame(
    {
        "Taxon": ancom_raw["Taxon"],
        "Method": "ANCOM_BC2",
        "Primary_Statistic": "Robust_BH_call",
        "Tested": ancom_tested,
        "P_Value": ancom_p,
        "Adjusted_P_Value": ancom_q,
        "Significant": ancom_primary_sig,
        "Effect_Estimate": ancom_effect,
        "Effect_Metric": "ANCOM-BC2 log fold change",
        "Effect_Direction": ancom_direction,
    }
)


# ============================================================
# 10. STANDARDIZE DESEQ2
# ============================================================

deseq_p_col = find_column(
    deseq_raw,
    [
        "P_Value",
        "pvalue",
    ]
)

deseq_q_col = find_column(
    deseq_raw,
    [
        "Adjusted_P_Value",
        "padj",
    ]
)

deseq_sig_col = find_column(
    deseq_raw,
    [
        "Primary_Significant",
        "Significant",
    ],
    required=False
)

deseq_tested_col = find_column(
    deseq_raw,
    [
        "Tested",
    ],
    required=False
)

deseq_effect_col = find_column(
    deseq_raw,
    [
        "Log2_Fold_Change",
        "log2FoldChange",
    ],
    required=False
)

deseq_direction_col = find_column(
    deseq_raw,
    [
        "Effect_Direction",
        "Direction",
    ],
    required=False
)


deseq_p = numeric_series(
    deseq_raw,
    deseq_p_col
)

deseq_q = numeric_series(
    deseq_raw,
    deseq_q_col
)


if deseq_sig_col is not None:

    deseq_supplied_sig = boolean_series(
        deseq_raw,
        deseq_sig_col
    )

else:

    deseq_supplied_sig = None


deseq_sig = derive_significance(
    deseq_q,
    deseq_supplied_sig
)


if deseq_tested_col is not None:

    deseq_tested = boolean_series(
        deseq_raw,
        deseq_tested_col
    )

else:

    deseq_tested = (
        deseq_p.notna()
    ).astype("boolean")


deseq_effect = numeric_series(
    deseq_raw,
    deseq_effect_col
)


if deseq_direction_col is not None:

    deseq_direction = (
        deseq_raw[deseq_direction_col]
        .astype("string")
    )

else:

    deseq_direction = pd.Series(
        pd.NA,
        index=deseq_raw.index,
        dtype="string"
    )


deseq_long = pd.DataFrame(
    {
        "Taxon": deseq_raw["Taxon"],
        "Method": "DESeq2",
        "Primary_Statistic": "Wald",
        "Tested": deseq_tested,
        "P_Value": deseq_p,
        "Adjusted_P_Value": deseq_q,
        "Significant": deseq_sig,
        "Effect_Estimate": deseq_effect,
        "Effect_Metric": "DESeq2 log2 fold change",
        "Effect_Direction": deseq_direction,
    }
)


# ============================================================
# 11. COMBINE FOUR METHODS INTO LONG FORMAT
# ============================================================

long_df = pd.concat(
    [
        clr_long,
        aldex_long,
        ancom_long,
        deseq_long,
    ],
    ignore_index=True
)


# ============================================================
# 12. ADD TRUTH INFORMATION
# ============================================================

truth_columns = [
    column
    for column in [
        "Taxon",
        "True_DA",
        "True_Direction",
        "True_Effect_Size",
    ]
    if column in truth_raw.columns
]

truth_for_merge = truth_raw[
    truth_columns
].copy()

long_df = long_df.merge(
    truth_for_merge,
    on="Taxon",
    how="left",
    validate="many_to_one"
)

long_df["Scenario"] = (
    "Null_Baseline_Low_Sparsity"
)

long_df["Seed"] = 42

long_df = long_df[
    [
        "Taxon",
        "Method",
        "Primary_Statistic",
        "Tested",
        "P_Value",
        "Adjusted_P_Value",
        "Significant",
        "Effect_Estimate",
        "Effect_Metric",
        "Effect_Direction",
    ]
    +
    [
        column
        for column in [
            "True_DA",
            "True_Direction",
            "True_Effect_Size",
        ]
        if column in long_df.columns
    ]
    +
    [
        "Scenario",
        "Seed",
    ]
]


# ============================================================
# 13. QC LONG TABLE
# ============================================================

expected_rows = (
    100 * 4
)

if len(long_df) != expected_rows:

    raise ValueError(
        "Expected 400 rows in standardized long table; "
        f"found {len(long_df)}."
    )

duplicates = long_df.duplicated(
    subset=[
        "Taxon",
        "Method",
    ]
).sum()

if duplicates != 0:

    raise ValueError(
        f"Found {duplicates} duplicated Taxon/Method rows."
    )


method_counts = (
    long_df.groupby("Method")["Taxon"]
    .nunique()
)

if not (
    method_counts == 100
).all():

    raise ValueError(
        "Not every method contains exactly 100 unique taxa:\n"
        f"{method_counts}"
    )


# ============================================================
# 14. METHOD SUMMARY
# ============================================================

method_summary_rows = []

for method, method_df in long_df.groupby(
    "Method",
    sort=False
):

    tested_count = int(
        method_df["Tested"]
        .fillna(False)
        .sum()
    )

    significant_count = int(
        method_df["Significant"]
        .fillna(False)
        .sum()
    )

    adjusted_available = int(
        method_df[
            "Adjusted_P_Value"
        ].notna().sum()
    )

    if significant_count > 0:
        realized_fdp = 1.0
    else:
        realized_fdp = 0.0

    primary_statistic = (
        method_df[
            "Primary_Statistic"
        ].iloc[0]
    )

    method_summary_rows.append(
        {
            "Method": method,
            "Primary_Statistic": (
                primary_statistic
            ),
            "Taxa_Available": (
                method_df["Taxon"].nunique()
            ),
            "Taxa_Tested": tested_count,
            "Adjusted_P_Available": (
                adjusted_available
            ),
            "Significant_Taxa": (
                significant_count
            ),
            "Minimum_Raw_P": safe_min(
                method_df["P_Value"]
            ),
            "Minimum_Adjusted_P": safe_min(
                method_df[
                    "Adjusted_P_Value"
                ]
            ),
            "Known_False_Positive_Discoveries": (
                significant_count
            ),
            "Realized_FDP_This_Null_Replicate": (
                realized_fdp
            ),
        }
    )

method_summary = pd.DataFrame(
    method_summary_rows
)


# ============================================================
# 15. CREATE WIDE FOUR-METHOD TABLE
# ============================================================

wide = truth_raw.copy()

wide = wide[
    [
        column
        for column in [
            "Taxon",
            "True_DA",
            "True_Direction",
            "True_Effect_Size",
        ]
        if column in wide.columns
    ]
].copy()


def add_method_to_wide(
    wide_df,
    method_df,
    prefix
):

    selected = method_df[
        [
            "Taxon",
            "Tested",
            "P_Value",
            "Adjusted_P_Value",
            "Significant",
            "Effect_Estimate",
            "Effect_Metric",
            "Effect_Direction",
        ]
    ].copy()

    selected = selected.rename(
        columns={
            "Tested": f"{prefix}_Tested",
            "P_Value": f"{prefix}_P",
            "Adjusted_P_Value": f"{prefix}_Q",
            "Significant": f"{prefix}_Significant",
            "Effect_Estimate": (
                f"{prefix}_Effect"
            ),
            "Effect_Metric": (
                f"{prefix}_Effect_Metric"
            ),
            "Effect_Direction": (
                f"{prefix}_Direction"
            ),
        }
    )

    return wide_df.merge(
        selected,
        on="Taxon",
        how="left",
        validate="one_to_one"
    )


wide = add_method_to_wide(
    wide,
    clr_long,
    "CLR"
)

wide = add_method_to_wide(
    wide,
    aldex_long,
    "ALDEx2"
)

wide = add_method_to_wide(
    wide,
    ancom_long,
    "ANCOM_BC2"
)

wide = add_method_to_wide(
    wide,
    deseq_long,
    "DESeq2"
)


# ============================================================
# 16. ADD ALDEX2 WELCH DIAGNOSTICS
# ============================================================

aldex_welch_p_col = find_column(
    aldex_raw,
    [
        "Welch_P_Value",
        "we.ep",
    ],
    required=False
)

aldex_welch_q_col = find_column(
    aldex_raw,
    [
        "Welch_Adjusted_P_Value",
        "we.eBH",
    ],
    required=False
)

aldex_welch_sig_col = find_column(
    aldex_raw,
    [
        "Welch_Significant",
    ],
    required=False
)

if (
    aldex_welch_p_col is not None and
    aldex_welch_q_col is not None
):

    welch_diag = pd.DataFrame(
        {
            "Taxon": aldex_raw["Taxon"],
            "ALDEx2_Welch_P": numeric_series(
                aldex_raw,
                aldex_welch_p_col
            ),
            "ALDEx2_Welch_Q": numeric_series(
                aldex_raw,
                aldex_welch_q_col
            ),
        }
    )

    if aldex_welch_sig_col is not None:

        welch_diag[
            "ALDEx2_Welch_Significant"
        ] = boolean_series(
            aldex_raw,
            aldex_welch_sig_col
        )

    else:

        welch_diag[
            "ALDEx2_Welch_Significant"
        ] = (
            welch_diag[
                "ALDEx2_Welch_Q"
            ] < ALPHA
        ).astype("boolean")

    wide = wide.merge(
        welch_diag,
        on="Taxon",
        how="left",
        validate="one_to_one"
    )


# ============================================================
# 17. ADD ANCOM STANDARD VS ROBUST DIAGNOSTICS
# ============================================================

ancom_diag = pd.DataFrame(
    {
        "Taxon": ancom_raw["Taxon"],
        "ANCOM_BC2_Standard_Significant": (
            ancom_standard_sig
        ),
        "ANCOM_BC2_Robust_Significant": (
            ancom_primary_sig
        ),
    }
)

wide = wide.merge(
    ancom_diag,
    on="Taxon",
    how="left",
    validate="one_to_one"
)


# ============================================================
# 18. CONSENSUS VOTING
# ============================================================

significance_columns = {
    "CLR_Wilcoxon": "CLR_Significant",
    "ALDEx2": "ALDEx2_Significant",
    "ANCOM_BC2": "ANCOM_BC2_Significant",
    "DESeq2": "DESeq2_Significant",
}

direction_columns = {
    "CLR_Wilcoxon": "CLR_Direction",
    "ALDEx2": "ALDEx2_Direction",
    "ANCOM_BC2": "ANCOM_BC2_Direction",
    "DESeq2": "DESeq2_Direction",
}


for column in significance_columns.values():

    wide[column] = (
        wide[column]
        .fillna(False)
        .astype(bool)
    )


wide["Significant_Method_Count"] = sum(
    wide[column].astype(int)
    for column
    in significance_columns.values()
)

wide["Consensus_2of4"] = (
    wide["Significant_Method_Count"] >= 2
)

wide["Consensus_3of4"] = (
    wide["Significant_Method_Count"] >= 3
)

wide["Consensus_4of4"] = (
    wide["Significant_Method_Count"] == 4
)


# ============================================================
# 19. DIRECTION-AWARE CONSENSUS
# ============================================================

def evaluate_direction_consensus(row):

    directions = []

    supporting_methods = []

    for method in significance_columns:

        sig_column = (
            significance_columns[method]
        )

        direction_column = (
            direction_columns[method]
        )

        if bool(row[sig_column]):

            supporting_methods.append(
                method
            )

            direction = row[
                direction_column
            ]

            if (
                pd.notna(direction) and
                str(direction).strip() not in {
                    "",
                    "<NA>",
                    "nan",
                    "None",
                }
            ):
                directions.append(
                    str(direction)
                )

    vote_count = len(
        supporting_methods
    )

    if vote_count == 0:

        return pd.Series(
            {
                "Supporting_Methods": "",
                "Significant_Direction_Agreement": (
                    pd.NA
                ),
                "Consensus_Direction": (
                    "No_significant_methods"
                ),
            }
        )

    unique_directions = sorted(
        set(directions)
    )

    if (
        len(directions) == vote_count and
        len(unique_directions) == 1
    ):

        agreement = True

        consensus_direction = (
            unique_directions[0]
        )

    elif len(unique_directions) > 1:

        agreement = False

        consensus_direction = (
            "Discordant"
        )

    else:

        agreement = pd.NA

        consensus_direction = (
            "Direction_incomplete"
        )

    return pd.Series(
        {
            "Supporting_Methods": ";".join(
                supporting_methods
            ),
            "Significant_Direction_Agreement": (
                agreement
            ),
            "Consensus_Direction": (
                consensus_direction
            ),
        }
    )


direction_result = wide.apply(
    evaluate_direction_consensus,
    axis=1
)

wide = pd.concat(
    [
        wide,
        direction_result,
    ],
    axis=1
)


wide["Direction_Aware_2of4"] = (
    wide["Consensus_2of4"] &
    (
        wide[
            "Significant_Direction_Agreement"
        ] == True
    )
)

wide["Direction_Aware_3of4"] = (
    wide["Consensus_3of4"] &
    (
        wide[
            "Significant_Direction_Agreement"
        ] == True
    )
)

wide["Direction_Aware_4of4"] = (
    wide["Consensus_4of4"] &
    (
        wide[
            "Significant_Direction_Agreement"
        ] == True
    )
)


# ============================================================
# 20. CREATE CONSENSUS TABLE
# ============================================================

consensus_columns = [
    "Taxon",
]

for column in [
    "True_DA",
    "True_Direction",
    "True_Effect_Size",
]:
    if column in wide.columns:
        consensus_columns.append(
            column
        )

consensus_columns.extend(
    [
        "CLR_Significant",
        "CLR_Direction",
        "ALDEx2_Significant",
        "ALDEx2_Direction",
        "ANCOM_BC2_Significant",
        "ANCOM_BC2_Direction",
        "DESeq2_Significant",
        "DESeq2_Direction",
        "Significant_Method_Count",
        "Supporting_Methods",
        "Consensus_2of4",
        "Consensus_3of4",
        "Consensus_4of4",
        "Significant_Direction_Agreement",
        "Consensus_Direction",
        "Direction_Aware_2of4",
        "Direction_Aware_3of4",
        "Direction_Aware_4of4",
    ]
)

consensus_df = wide[
    consensus_columns
].copy()


# ============================================================
# 21. PAIRWISE JACCARD OF SIGNIFICANT SETS
#
# If both methods have zero discoveries, Jaccard is undefined
# and is reported as NA rather than artificially assigning 1.
# ============================================================

method_names = list(
    significance_columns.keys()
)

significant_sets = {}

for method in method_names:

    column = significance_columns[
        method
    ]

    significant_sets[method] = set(
        wide.loc[
            wide[column],
            "Taxon"
        ]
    )


jaccard = pd.DataFrame(
    index=method_names,
    columns=method_names,
    dtype=float
)

for method_a in method_names:

    for method_b in method_names:

        if method_a == method_b:

            jaccard.loc[
                method_a,
                method_b
            ] = 1.0

            continue

        set_a = significant_sets[
            method_a
        ]

        set_b = significant_sets[
            method_b
        ]

        union = (
            set_a | set_b
        )

        intersection = (
            set_a & set_b
        )

        if len(union) == 0:

            value = np.nan

        else:

            value = (
                len(intersection) /
                len(union)
            )

        jaccard.loc[
            method_a,
            method_b
        ] = value


# ============================================================
# 22. PAIRWISE SIGNIFICANCE-CALL AGREEMENT
#
# This measures agreement across all 100 taxa, not just
# overlap among positive discoveries.
# ============================================================

call_agreement = pd.DataFrame(
    index=method_names,
    columns=method_names,
    dtype=float
)

for method_a in method_names:

    for method_b in method_names:

        column_a = significance_columns[
            method_a
        ]

        column_b = significance_columns[
            method_b
        ]

        agreement = (
            wide[column_a] ==
            wide[column_b]
        ).mean()

        call_agreement.loc[
            method_a,
            method_b
        ] = float(
            agreement
        )


# ============================================================
# 23. TAXON_033 FOCUSED TABLE
# ============================================================

taxon033_long = long_df[
    long_df["Taxon"] ==
    "Taxon_033"
].copy()

if len(taxon033_long) != 4:

    raise ValueError(
        "Expected four Taxon_033 method rows; "
        f"found {len(taxon033_long)}."
    )


# ============================================================
# 24. FINAL QC
# ============================================================

total_2of4 = int(
    wide["Consensus_2of4"].sum()
)

total_3of4 = int(
    wide["Consensus_3of4"].sum()
)

total_4of4 = int(
    wide["Consensus_4of4"].sum()
)

total_direction_2of4 = int(
    wide[
        "Direction_Aware_2of4"
    ].sum()
)


print("\nMethod summary:")

print(
    method_summary.to_string(
        index=False
    )
)


print("\nConsensus summary:")

print(
    f"Taxa supported by >=2 methods: "
    f"{total_2of4}"
)

print(
    f"Taxa supported by >=3 methods: "
    f"{total_3of4}"
)

print(
    f"Taxa supported by 4 methods: "
    f"{total_4of4}"
)

print(
    f"Direction-aware >=2-method taxa: "
    f"{total_direction_2of4}"
)


print("\nTaxon_033 comparison:")

print(
    taxon033_long[
        [
            "Taxon",
            "Method",
            "Primary_Statistic",
            "P_Value",
            "Adjusted_P_Value",
            "Significant",
            "Effect_Estimate",
            "Effect_Direction",
        ]
    ].to_string(
        index=False
    )
)


# ============================================================
# 25. SAVE OUTPUTS
# ============================================================

long_df.to_csv(
    LONG_FILE,
    index=False
)

wide.to_csv(
    WIDE_FILE,
    index=False
)

method_summary.to_csv(
    METHOD_SUMMARY_FILE,
    index=False
)

consensus_df.to_csv(
    CONSENSUS_FILE,
    index=False
)

taxon033_long.to_csv(
    TAXON033_FILE,
    index=False
)

jaccard.to_csv(
    JACCARD_FILE,
    index=True
)

call_agreement.to_csv(
    CALL_AGREEMENT_FILE,
    index=True
)


# ============================================================
# 26. WRITE HUMAN-READABLE SUMMARY
# ============================================================

summary_lines = [
    "CONSOLIDATED FOUR-METHOD NULL VALIDATION",
    "=" * 64,
    "",
    "Dataset:",
    "100 taxa x 50 samples",
    "25 Control / 25 Case",
    "Known true DA taxa: 0",
    "",
    "Primary method conventions:",
    "CLR/Wilcoxon: BH-adjusted Wilcoxon result",
    (
        "ALDEx2: Wilcoxon expected BH result; "
        "Welch retained diagnostically"
    ),
    (
        "ANCOM-BC2: robust significance call "
        "when available"
    ),
    (
        "DESeq2: primary Wald result with "
        "standard independent filtering"
    ),
    "",
    "METHOD SUMMARY",
]

for _, row in method_summary.iterrows():

    summary_lines.append(
        (
            f"{row['Method']}: "
            f"tested={int(row['Taxa_Tested'])}, "
            f"significant={int(row['Significant_Taxa'])}, "
            f"min raw p={row['Minimum_Raw_P']:.10g}, "
            f"min adjusted p="
            f"{row['Minimum_Adjusted_P']:.10g}"
        )
    )


summary_lines.extend(
    [
        "",
        "CONSENSUS",
        (
            f">=2-of-4 significant: "
            f"{total_2of4}"
        ),
        (
            f">=3-of-4 significant: "
            f"{total_3of4}"
        ),
        (
            f"4-of-4 significant: "
            f"{total_4of4}"
        ),
        (
            "Direction-aware >=2-of-4: "
            f"{total_direction_2of4}"
        ),
        "",
        "INTERPRETATION",
        (
            "Because this simulation contains zero true "
            "differentially abundant taxa, every discovery "
            "is a known false-positive discovery."
        ),
        (
            "Agreement among multiple methods in this null "
            "replicate therefore does not establish truth."
        ),
        (
            "Single-replicate results describe this realization "
            "only; calibration must be estimated across repeated "
            "simulations."
        ),
        "",
        "OUTPUT FILES",
        str(LONG_FILE),
        str(WIDE_FILE),
        str(METHOD_SUMMARY_FILE),
        str(CONSENSUS_FILE),
        str(TAXON033_FILE),
        str(JACCARD_FILE),
        str(CALL_AGREEMENT_FILE),
    ]
)

SUMMARY_FILE.write_text(
    "\n".join(summary_lines),
    encoding="utf-8"
)


# ============================================================
# 27. COMPLETE
# ============================================================

print("\n" + "=" * 72)

print(
    "CONSOLIDATED NULL VALIDATION: PASS"
)

print("=" * 72)

print("\nOutput files:")

for output_file in [
    LONG_FILE,
    WIDE_FILE,
    METHOD_SUMMARY_FILE,
    CONSENSUS_FILE,
    TAXON033_FILE,
    JACCARD_FILE,
    CALL_AGREEMENT_FILE,
    SUMMARY_FILE,
]:
    print(output_file)
