import numpy as np
import pandas as pd


# ============================================================
# INPUT FILES
# ============================================================

COUNT_FILE = "data/processed/simulation_null/null_counts.csv"
METADATA_FILE = "data/processed/simulation_null/null_metadata.csv"
TRUTH_FILE = "data/processed/simulation_null/null_ground_truth.csv"


# ============================================================
# LOAD DATA
# ============================================================

counts = pd.read_csv(
    COUNT_FILE,
    index_col="Taxon"
)

metadata = pd.read_csv(
    METADATA_FILE
)

# IMPORTANT:
# keep_default_na=False prevents the literal string "None"
# from being interpreted as a missing value (NaN).
truth = pd.read_csv(
    TRUTH_FILE,
    keep_default_na=False
)


# ============================================================
# 1. DIMENSIONS
# ============================================================

print("=" * 70)
print("1. DIMENSIONS")
print("=" * 70)

print("Count table shape:", counts.shape)
print("Metadata shape:", metadata.shape)
print("Truth table shape:", truth.shape)


# ============================================================
# 2. IDENTIFIERS
# ============================================================

print()
print("=" * 70)
print("2. IDENTIFIERS")
print("=" * 70)

print(
    "Duplicate taxa:",
    counts.index.duplicated().sum()
)

print(
    "Duplicate metadata SampleIDs:",
    metadata["SampleID"].duplicated().sum()
)

print(
    "Duplicate truth taxa:",
    truth["Taxon"].duplicated().sum()
)

print()

print(
    "Count sample IDs match metadata IDs:",
    list(counts.columns) == list(metadata["SampleID"])
)

print(
    "Count taxa match truth taxa:",
    list(counts.index) == list(truth["Taxon"])
)


# ============================================================
# 3. GROUP STRUCTURE
# ============================================================

print()
print("=" * 70)
print("3. GROUP STRUCTURE")
print("=" * 70)

print(metadata["Group"].value_counts())


# ============================================================
# 4. COUNT INTEGRITY
# ============================================================

print()
print("=" * 70)
print("4. COUNT INTEGRITY")
print("=" * 70)

values = counts.to_numpy()

all_numeric = all(
    pd.api.types.is_numeric_dtype(dtype)
    for dtype in counts.dtypes
)

negative_counts = int(
    (values < 0).sum()
)

non_integer_values = int(
    np.sum(values != np.floor(values))
)

print(
    "All numeric:",
    all_numeric
)

print(
    "Negative counts:",
    negative_counts
)

print(
    "Non-integer values:",
    non_integer_values
)


# ============================================================
# 5. LIBRARY SIZE VALIDATION
# ============================================================

print()
print("=" * 70)
print("5. LIBRARY SIZE VALIDATION")
print("=" * 70)

observed_library_sizes = counts.sum(axis=0)

stored_library_sizes = (
    metadata
    .set_index("SampleID")["LibrarySize"]
    .reindex(counts.columns)
)

library_sizes_match = (
    observed_library_sizes
    == stored_library_sizes
).all()

print(
    "Library sizes match metadata:",
    bool(library_sizes_match)
)

print(
    "Mean observed library size:",
    round(observed_library_sizes.mean(), 2)
)

print(
    "Minimum observed library size:",
    int(observed_library_sizes.min())
)

print(
    "Maximum observed library size:",
    int(observed_library_sizes.max())
)

print(
    "Zero-library samples:",
    int((observed_library_sizes == 0).sum())
)


# ============================================================
# 6. SPARSITY
# ============================================================

print()
print("=" * 70)
print("6. SPARSITY")
print("=" * 70)

zero_cells = int(
    (values == 0).sum()
)

total_cells = values.size

zero_fraction = (
    zero_cells
    / total_cells
)

print(
    "Total cells:",
    total_cells
)

print(
    "Zero cells:",
    zero_cells
)

print(
    "Zero fraction:",
    round(zero_fraction, 4)
)


# ============================================================
# 7. TAXON PREVALENCE
# ============================================================

print()
print("=" * 70)
print("7. TAXON PREVALENCE")
print("=" * 70)

prevalence_counts = (
    counts > 0
).sum(axis=1)

prevalence_fraction = (
    prevalence_counts
    / counts.shape[1]
)

total_taxon_counts = counts.sum(axis=1)

print(
    "Minimum taxon prevalence:",
    int(prevalence_counts.min()),
    "samples"
)

print(
    "Median taxon prevalence:",
    float(prevalence_counts.median()),
    "samples"
)

print(
    "Maximum taxon prevalence:",
    int(prevalence_counts.max()),
    "samples"
)

all_zero_taxa = (
    prevalence_counts[
        prevalence_counts == 0
    ]
    .index
    .tolist()
)

print(
    "All-zero taxa:",
    len(all_zero_taxa)
)

print(
    "All-zero taxon IDs:",
    all_zero_taxa
)


# ============================================================
# 8. RAREST TAXA
# ============================================================

print()
print("=" * 70)
print("8. RAREST TAXA")
print("=" * 70)

prevalence_table = pd.DataFrame(
    {
        "Present_Samples": prevalence_counts,
        "Prevalence_Fraction": prevalence_fraction,
        "Total_Count": total_taxon_counts
    }
)

prevalence_table = prevalence_table.sort_values(
    by=[
        "Present_Samples",
        "Total_Count"
    ]
)

print(
    prevalence_table
    .head(15)
    .to_string()
)


# ============================================================
# 9. GROUND TRUTH
# ============================================================

print()
print("=" * 70)
print("9. GROUND TRUTH")
print("=" * 70)

print(
    "True DA taxa:",
    int(truth["True_DA"].sum())
)

print(
    "Unique true directions:",
    truth["True_Direction"].unique()
)

print(
    "Unique true effect sizes:",
    truth["True_Effect_Size"].unique()
)


# ============================================================
# 10. SUMMARY
# ============================================================

print()
print("=" * 70)
print("10. SUMMARY")
print("=" * 70)

print(
    "Mean count per cell:",
    round(values.mean(), 2)
)

print(
    "Median count per cell:",
    float(np.median(values))
)

print()

if (
    counts.shape == (100, 50)
    and metadata.shape[0] == 50
    and truth.shape[0] == 100
    and counts.index.duplicated().sum() == 0
    and metadata["SampleID"].duplicated().sum() == 0
    and truth["Taxon"].duplicated().sum() == 0
    and list(counts.columns) == list(metadata["SampleID"])
    and list(counts.index) == list(truth["Taxon"])
    and negative_counts == 0
    and non_integer_values == 0
    and bool(library_sizes_match)
    and int(truth["True_DA"].sum()) == 0
):
    print(
        "OVERALL RESULT: PASS — "
        "all critical integrity checks passed."
    )
else:
    print(
        "OVERALL RESULT: REVIEW REQUIRED — "
        "one or more critical checks failed."
    )

print()
print("Manual audit completed.")
