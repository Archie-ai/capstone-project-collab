from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# PATH CONFIGURATION
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data" / "processed" / "simulation_null"

COUNT_FILE = DATA_DIR / "null_counts.csv"
METADATA_FILE = DATA_DIR / "null_metadata.csv"
TRUTH_FILE = DATA_DIR / "null_ground_truth.csv"

QC_SUMMARY_FILE = DATA_DIR / "null_qc_summary.txt"
TAXON_QC_FILE = DATA_DIR / "null_taxon_qc.csv"


# ============================================================
# EXPECTED DESIGN
# ============================================================

EXPECTED_SAMPLES = 50
EXPECTED_TAXA = 100
EXPECTED_CONTROL = 25
EXPECTED_CASE = 25
EXPECTED_TRUE_DA = 0


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def report_check(name, passed, details=""):
    status = "PASS" if passed else "FAIL"

    message = f"[{status}] {name}"

    if details:
        message += f": {details}"

    print(message)

    return message, passed


def report_warning(name, condition, details=""):
    if condition:
        message = f"[WARNING] {name}"

        if details:
            message += f": {details}"

        print(message)

        return message

    return None


# ============================================================
# LOAD DATA
# ============================================================

print()
print("NULL DATASET QUALITY-CONTROL AUDIT")
print("=" * 50)
print()

count_df = pd.read_csv(
    COUNT_FILE,
    index_col="Taxon"
)

metadata_df = pd.read_csv(
    METADATA_FILE
)

truth_df = pd.read_csv(
    TRUTH_FILE
)

audit_messages = []
critical_checks = []


# ============================================================
# 1. DATASET DIMENSIONS
# ============================================================

message, passed = report_check(
    "Count-table dimensions",
    count_df.shape == (EXPECTED_TAXA, EXPECTED_SAMPLES),
    f"observed={count_df.shape}, expected=(100, 50)"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 2. DUPLICATE TAXA
# ============================================================

duplicate_taxa = count_df.index.duplicated().sum()

message, passed = report_check(
    "Unique taxon identifiers",
    duplicate_taxa == 0,
    f"duplicate taxa={duplicate_taxa}"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 3. DUPLICATE SAMPLE COLUMNS
# ============================================================

duplicate_sample_columns = count_df.columns.duplicated().sum()

message, passed = report_check(
    "Unique count-table sample IDs",
    duplicate_sample_columns == 0,
    f"duplicate sample columns={duplicate_sample_columns}"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 4. COUNT VALUES ARE NUMERIC
# ============================================================

numeric_counts = count_df.apply(
    lambda column: pd.to_numeric(
        column,
        errors="coerce"
    )
)

missing_numeric = numeric_counts.isna().sum().sum()

message, passed = report_check(
    "All count values are numeric",
    missing_numeric == 0,
    f"non-numeric/missing values={missing_numeric}"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 5. COUNTS ARE NON-NEGATIVE
# ============================================================

negative_counts = (numeric_counts < 0).sum().sum()

message, passed = report_check(
    "Counts are non-negative",
    negative_counts == 0,
    f"negative cells={negative_counts}"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 6. COUNTS ARE INTEGERS
# ============================================================

count_values = numeric_counts.to_numpy()

integer_counts = np.all(
    np.equal(
        count_values,
        np.floor(count_values)
    )
)

message, passed = report_check(
    "Counts are integers",
    integer_counts
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 7. METADATA SAMPLE COUNT
# ============================================================

message, passed = report_check(
    "Metadata sample count",
    len(metadata_df) == EXPECTED_SAMPLES,
    f"observed={len(metadata_df)}, expected=50"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 8. UNIQUE METADATA SAMPLE IDs
# ============================================================

duplicate_metadata_samples = metadata_df["SampleID"].duplicated().sum()

message, passed = report_check(
    "Unique metadata sample IDs",
    duplicate_metadata_samples == 0,
    f"duplicates={duplicate_metadata_samples}"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 9. COUNT / METADATA SAMPLE ALIGNMENT
# ============================================================

count_samples = list(count_df.columns)
metadata_samples = list(metadata_df["SampleID"])

sample_sets_match = set(count_samples) == set(metadata_samples)
sample_order_matches = count_samples == metadata_samples

message, passed = report_check(
    "Count-table and metadata sample IDs match",
    sample_sets_match
)

audit_messages.append(message)
critical_checks.append(passed)

message, passed = report_check(
    "Count-table and metadata sample order match",
    sample_order_matches
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 10. GROUP BALANCE
# ============================================================

group_counts = metadata_df["Group"].value_counts()

control_count = int(group_counts.get("Control", 0))
case_count = int(group_counts.get("Case", 0))

group_balance_ok = (
    control_count == EXPECTED_CONTROL
    and case_count == EXPECTED_CASE
)

message, passed = report_check(
    "Case/Control group balance",
    group_balance_ok,
    f"Control={control_count}, Case={case_count}"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 11. LIBRARY SIZE VALIDATION
# ============================================================

observed_library_sizes = count_df.sum(axis=0)

metadata_library_sizes = metadata_df.set_index(
    "SampleID"
)["LibrarySize"]

metadata_library_sizes = metadata_library_sizes.loc[
    count_df.columns
]

library_size_match = np.array_equal(
    observed_library_sizes.to_numpy(),
    metadata_library_sizes.to_numpy()
)

message, passed = report_check(
    "Recorded library sizes equal observed count sums",
    library_size_match
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 12. ZERO-LIBRARY SAMPLES
# ============================================================

zero_library_samples = int(
    (observed_library_sizes == 0).sum()
)

message, passed = report_check(
    "No zero-library samples",
    zero_library_samples == 0,
    f"zero-library samples={zero_library_samples}"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 13. ZERO FRACTION
# ============================================================

zero_fraction = (
    (count_df == 0).sum().sum()
    / count_df.size
)

print(
    f"[INFO] Overall zero fraction: "
    f"{zero_fraction:.4f} "
    f"({zero_fraction * 100:.2f}%)"
)

audit_messages.append(
    f"[INFO] Overall zero fraction: "
    f"{zero_fraction:.4f} "
    f"({zero_fraction * 100:.2f}%)"
)


# ============================================================
# 14. TAXON PREVALENCE
# ============================================================

taxon_prevalence = (
    (count_df > 0).sum(axis=1)
    / count_df.shape[1]
)

print(
    "[INFO] Taxon prevalence: "
    f"min={taxon_prevalence.min():.3f}, "
    f"median={taxon_prevalence.median():.3f}, "
    f"max={taxon_prevalence.max():.3f}"
)

audit_messages.append(
    "[INFO] Taxon prevalence: "
    f"min={taxon_prevalence.min():.3f}, "
    f"median={taxon_prevalence.median():.3f}, "
    f"max={taxon_prevalence.max():.3f}"
)


# ============================================================
# 15. ALL-ZERO TAXA
# ============================================================

all_zero_mask = count_df.sum(axis=1) == 0
all_zero_taxa = count_df.index[all_zero_mask].tolist()

warning = report_warning(
    "All-zero taxa detected",
    len(all_zero_taxa) > 0,
    (
        f"count={len(all_zero_taxa)}, "
        f"taxa={', '.join(all_zero_taxa)}"
    )
)

if warning:
    audit_messages.append(warning)

print(
    f"[INFO] Number of all-zero taxa: "
    f"{len(all_zero_taxa)}"
)

audit_messages.append(
    f"[INFO] Number of all-zero taxa: "
    f"{len(all_zero_taxa)}"
)


# ============================================================
# 16. GROUND-TRUTH DIMENSIONS
# ============================================================

message, passed = report_check(
    "Ground-truth taxon count",
    len(truth_df) == EXPECTED_TAXA,
    f"observed={len(truth_df)}, expected=100"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 17. TAXON / TRUTH ALIGNMENT
# ============================================================

truth_taxa = list(truth_df["Taxon"])

truth_alignment = set(count_df.index) == set(truth_taxa)

message, passed = report_check(
    "Count-table and ground-truth taxa match",
    truth_alignment
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 18. TRUE DA STATUS
# ============================================================

true_da_count = int(
    truth_df["True_DA"].sum()
)

message, passed = report_check(
    "Null ground truth contains zero DA taxa",
    true_da_count == EXPECTED_TRUE_DA,
    f"true DA taxa={true_da_count}"
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 19. TRUE EFFECT SIZE
# ============================================================

effect_size_zero = bool(
    (truth_df["True_Effect_Size"] == 0).all()
)

message, passed = report_check(
    "All null effect sizes equal zero",
    effect_size_zero
)

audit_messages.append(message)
critical_checks.append(passed)


# ============================================================
# 20. GROUP LIBRARY-SIZE SUMMARY
# ============================================================

library_summary_df = metadata_df.groupby(
    "Group"
)["LibrarySize"].agg(
    ["mean", "median", "min", "max"]
)

print()
print("Library-size summary by group:")
print(library_summary_df.round(2))

audit_messages.append("")
audit_messages.append(
    "Library-size summary by group:"
)

audit_messages.append(
    library_summary_df.round(2).to_string()
)


# ============================================================
# 21. CREATE TAXON-LEVEL QC TABLE
# ============================================================

taxon_qc_df = pd.DataFrame(
    {
        "Taxon": count_df.index,
        "Total_Count": count_df.sum(axis=1).values,
        "Mean_Count": count_df.mean(axis=1).values,
        "Zero_Count": (count_df == 0).sum(axis=1).values,
        "Prevalence": taxon_prevalence.values,
        "All_Zero": all_zero_mask.values,
    }
)

taxon_qc_df.to_csv(
    TAXON_QC_FILE,
    index=False
)


# ============================================================
# FINAL AUDIT RESULT
# ============================================================

overall_pass = all(critical_checks)

print()
print("=" * 50)

if overall_pass:
    final_message = (
        "OVERALL RESULT: PASS — "
        "all critical integrity checks passed."
    )
else:
    final_message = (
        "OVERALL RESULT: FAIL — "
        "one or more critical checks failed."
    )

print(final_message)

audit_messages.append("")
audit_messages.append(final_message)

if all_zero_taxa:
    interpretation = (
        "Interpretation: the dataset passes structural QC, "
        "but all-zero taxa should be documented because "
        "different DA methods may filter or handle them differently."
    )

    print(interpretation)
    audit_messages.append(interpretation)


# ============================================================
# SAVE QC SUMMARY
# ============================================================

QC_SUMMARY_FILE.write_text(
    "\n".join(audit_messages) + "\n"
)

print()
print(f"QC summary saved to: {QC_SUMMARY_FILE}")
print(f"Taxon QC table saved to: {TAXON_QC_FILE}")

