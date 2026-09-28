#!/usr/bin/env python3

"""
Shared quality-control step for the microbiome consensus pipeline.

The goal is to give all four methods the same validated input.

- Counts stay as raw integers.
- No normalization happens here.
- No CLR transformation happens here.
- No pseudocount is added here.
- Ground truth is not used here.

Method-specific processing happens later inside each method.

This script was developed with the assitance from OpenAI's ChatGPT for code drafting, review, and documentation.
"""

import argparse
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd


# ============================================================
# COMMAND-LINE INPUTS
# ============================================================

def parse_args():
    """Collect the input files and QC settings from the command line."""

    parser = argparse.ArgumentParser(
        description="Shared QC before the four differential-abundance methods."
    )

    parser.add_argument(
        "--counts",
        required=True,
        type=Path,
        help="Raw count table: taxa in rows, samples in columns."
    )

    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Metadata file containing SampleID and Group."
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Folder where QC outputs will be saved."
    )

    parser.add_argument(
        "--dataset-label",
        required=True,
        help="Short dataset name used in output filenames, e.g. null or real_test."
    )

    parser.add_argument(
        "--sample-id-column",
        default="SampleID",
        help="Metadata column containing sample IDs. Default: SampleID."
    )

    parser.add_argument(
        "--group-column",
        default="Group",
        help="Metadata column containing case/control groups. Default: Group."
    )

    parser.add_argument(
        "--min-prevalence",
        type=float,
        default=0.0,
        help=(
            "Minimum fraction of samples where a taxon must be present. "
            "Default 0.0 means no prevalence filtering."
        )
    )

    parser.add_argument(
        "--min-total-count",
        type=int,
        default=1,
        help=(
            "Minimum total abundance across all samples. "
            "Default 1 removes only all-zero taxa."
        )
    )

    args = parser.parse_args()

    if not 0 <= args.min_prevalence <= 1:
        parser.error("--min-prevalence must be between 0 and 1.")

    if args.min_total_count < 1:
        parser.error("--min-total-count must be at least 1.")

    return args


# ============================================================
# FILE LOADING
# ============================================================

def load_data(count_file, metadata_file):
    """Load the count table and metadata."""

    if not count_file.exists():
        raise FileNotFoundError(f"Count file not found: {count_file}")

    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    # First column contains the taxon names.
    counts = pd.read_csv(count_file, index_col=0)
    metadata = pd.read_csv(metadata_file)

    if counts.empty:
        raise ValueError("Count table is empty.")

    if metadata.empty:
        raise ValueError("Metadata file is empty.")

    return counts, metadata


# ============================================================
# CHECK IDS AND METADATA
# ============================================================

def validate_ids(
    counts,
    metadata,
    sample_id_column,
    group_column
):
    """
    Check that taxa and samples have usable IDs and that the
    metadata contains the columns needed by the downstream methods.
    """

    required_columns = {sample_id_column, group_column}
    missing_columns = required_columns - set(metadata.columns)

    if missing_columns:
        raise ValueError(
            "Metadata is missing required column(s): "
            + ", ".join(sorted(missing_columns))
        )

    # Standardize IDs as strings and remove accidental spaces.
    counts = counts.copy()
    metadata = metadata.copy()

    counts.index = counts.index.astype(str).str.strip()
    counts.columns = counts.columns.astype(str).str.strip()

    metadata[sample_id_column] = (
        metadata[sample_id_column]
        .astype(str)
        .str.strip()
    )

    metadata[group_column] = (
        metadata[group_column]
        .astype(str)
        .str.strip()
    )

    if (counts.index == "").any():
        raise ValueError("Blank taxon ID found.")

    if (counts.columns == "").any():
        raise ValueError("Blank sample ID found in count table.")

    if (metadata[sample_id_column] == "").any():
        raise ValueError("Blank sample ID found in metadata.")

    if (metadata[group_column] == "").any():
        raise ValueError("Blank group label found in metadata.")

    if counts.index.duplicated().any():
        raise ValueError("Duplicate taxon IDs found.")

    if counts.columns.duplicated().any():
        raise ValueError("Duplicate sample IDs found in count table.")

    if metadata[sample_id_column].duplicated().any():
        raise ValueError("Duplicate sample IDs found in metadata.")

    if metadata[group_column].nunique() < 2:
        raise ValueError(
            "At least two groups are required for differential-abundance testing."
        )

    return counts, metadata


# ============================================================
# MATCH COUNTS TO METADATA
# ============================================================

def align_samples(counts, metadata, sample_id_column):
    """
    Make sure the same samples exist in both files.

    Metadata is reordered to match the count-table column order.
    This prevents a sample from accidentally receiving the wrong
    case/control label.
    """

    count_samples = list(counts.columns)
    metadata_samples = list(metadata[sample_id_column])

    only_in_counts = sorted(
        set(count_samples) - set(metadata_samples)
    )

    only_in_metadata = sorted(
        set(metadata_samples) - set(count_samples)
    )

    if only_in_counts or only_in_metadata:

        message = ["Sample IDs do not match between counts and metadata."]

        if only_in_counts:
            message.append(
                "Only in count table: "
                + ", ".join(only_in_counts[:10])
            )

        if only_in_metadata:
            message.append(
                "Only in metadata: "
                + ", ".join(only_in_metadata[:10])
            )

        raise ValueError("\n".join(message))

    # Put metadata rows in exactly the same order as count-table columns.
    metadata = (
        metadata
        .set_index(sample_id_column)
        .loc[count_samples]
        .reset_index()
    )

    return counts, metadata


# ============================================================
# CHECK RAW COUNTS
# ============================================================

def validate_counts(counts):
    """
    Confirm that the table contains valid raw sequencing counts.

    The shared QC output must remain integer counts because all four
    methods branch from this same dataset.
    """

    numeric_counts = counts.apply(
        pd.to_numeric,
        errors="coerce"
    )

    if numeric_counts.isna().any().any():
        raise ValueError(
            "Count table contains missing or non-numeric values."
        )

    values = numeric_counts.to_numpy(dtype=float)

    if not np.isfinite(values).all():
        raise ValueError("Count table contains infinite values.")

    if (values < 0).any():
        raise ValueError("Count table contains negative values.")

    # Values such as 5.0 are acceptable, but 5.3 is not a raw count.
    if not np.isclose(values, np.round(values)).all():
        raise ValueError(
            "Count table contains non-integer values."
        )

    return numeric_counts.round().astype(np.int64)


# ============================================================
# SAMPLE QC
# ============================================================

def make_sample_summary(
    counts,
    metadata,
    sample_id_column,
    group_column
):
    """Calculate basic QC information for each sample."""

    summary = pd.DataFrame({
        "SampleID": counts.columns,

        # Total sequencing depth for each sample.
        "LibrarySize": counts.sum(axis=0).values,

        # Number of taxa detected at least once.
        "DetectedTaxa": (counts > 0).sum(axis=0).values,

        # Number and fraction of taxa with zero counts.
        "ZeroCountTaxa": (counts == 0).sum(axis=0).values,
        "ZeroFraction": (counts == 0).mean(axis=0).values
    })

    group_lookup = (
        metadata
        .set_index(sample_id_column)[group_column]
    )

    summary[group_column] = (
        summary["SampleID"]
        .map(group_lookup)
    )

    # A completely empty sample cannot be used by the DA methods.
    summary["QC_Flag"] = np.where(
        summary["LibrarySize"] == 0,
        "FAIL_zero_library_size",
        "PASS"
    )

    return summary


# ============================================================
# TAXON QC
# ============================================================

def make_taxon_summary(
    counts,
    min_prevalence,
    min_total_count
):
    """Calculate QC information for every taxon."""

    total_count = counts.sum(axis=1)

    samples_detected = (
        counts > 0
    ).sum(axis=1)

    prevalence = (
        counts > 0
    ).mean(axis=1)

    zero_fraction = (
        counts == 0
    ).mean(axis=1)

    # A taxon must pass both selected thresholds.
    keep = (
        (total_count >= min_total_count)
        & (prevalence >= min_prevalence)
    )

    summary = pd.DataFrame({
        "Taxon": counts.index,
        "TotalCount": total_count.values,
        "SamplesDetected": samples_detected.values,
        "Prevalence": prevalence.values,
        "ZeroFraction": zero_fraction.values,
        "Keep": keep.values
    })

    # Make the reason for removal visible in the QC report.
    def qc_status(row):

        reasons = []

        if row["TotalCount"] < min_total_count:
            reasons.append("total_count_below_threshold")

        if row["Prevalence"] < min_prevalence:
            reasons.append("prevalence_below_threshold")

        if reasons:
            return ";".join(reasons)

        return "kept"

    summary["QC_Status"] = summary.apply(
        qc_status,
        axis=1
    )

    return summary


# ============================================================
# FILTER TAXA
# ============================================================

def filter_taxa(counts, taxon_summary):
    """
    Apply the shared taxon QC filter.

    With the default settings, this removes only taxa that are
    zero across every sample.
    """

    taxa_to_keep = taxon_summary.loc[
        taxon_summary["Keep"],
        "Taxon"
    ]

    filtered_counts = counts.loc[taxa_to_keep]

    if filtered_counts.empty:
        raise ValueError(
            "QC removed every taxon. Check the data or filtering thresholds."
        )

    return filtered_counts


# ============================================================
# WRITE HUMAN-READABLE SUMMARY
# ============================================================

def write_summary(
    output_file,
    dataset_label,
    original_counts,
    filtered_counts,
    metadata,
    sample_summary,
    taxon_summary,
    group_column,
    min_prevalence,
    min_total_count
):
    """Create a short text summary of the QC results."""

    taxa_before = original_counts.shape[0]
    taxa_after = filtered_counts.shape[0]

    zero_before = (
        original_counts == 0
    ).to_numpy().mean()

    zero_after = (
        filtered_counts == 0
    ).to_numpy().mean()

    group_counts = metadata[group_column].value_counts()

    group_text = "\n".join(
        f"  {group}: {count}"
        for group, count in group_counts.items()
    )

    failed_samples = (
        sample_summary["QC_Flag"] != "PASS"
    ).sum()

    all_zero_taxa = (
        taxon_summary["TotalCount"] == 0
    ).sum()

    report = f"""
MICROBIOME QUALITY-CONTROL SUMMARY
==================================

Dataset:
{dataset_label}

Samples:
{original_counts.shape[1]}

Groups:
{group_text}

Taxa before QC:
{taxa_before}

Taxa after QC:
{taxa_after}

Taxa removed:
{taxa_before - taxa_after}

All-zero taxa:
{all_zero_taxa}

Samples failing QC:
{failed_samples}

Mean library size:
{sample_summary["LibrarySize"].mean():.2f}

Median library size:
{sample_summary["LibrarySize"].median():.2f}

Minimum library size:
{sample_summary["LibrarySize"].min()}

Maximum library size:
{sample_summary["LibrarySize"].max()}

Zero fraction before QC:
{zero_before:.4f}

Zero fraction after QC:
{zero_after:.4f}

Minimum total-count threshold:
{min_total_count}

Minimum prevalence threshold:
{min_prevalence}

DOWNSTREAM INPUT
----------------
The QC count table contains raw integer counts.

The same QC counts and metadata are used as the starting input for:

1. CLR/Wilcoxon
2. ALDEx2
3. ANCOM-BC
4. DESeq2

No normalization, pseudocount addition, CLR transformation,
or differential-abundance testing was performed during QC.
"""

    output_file.write_text(
        report.strip() + "\n"
    )


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():

    args = parse_args()

    # Make labels safe for filenames:
    # "Real Test" becomes "real_test".
    dataset_label = (
        args.dataset_label
        .strip()
        .lower()
    )

    dataset_label = re.sub(
        r"[^a-z0-9]+",
        "_",
        dataset_label
    ).strip("_")

    if not dataset_label:
        raise ValueError(
            "Dataset label must contain letters or numbers."
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    print()
    print(f"Running QC for: {dataset_label}")
    print()

    # --------------------------------------------------------
    # 1. Load the input count table and metadata
    # --------------------------------------------------------

    counts, metadata = load_data(
        args.counts,
        args.metadata
    )

    # --------------------------------------------------------
    # 2. Check IDs and required metadata.
    # --------------------------------------------------------

    counts, metadata = validate_ids(
        counts,
        metadata,
        args.sample_id_column,
        args.group_column
    )

    # --------------------------------------------------------
    # 3. Match metadata to the count-table sample order.
    # --------------------------------------------------------

    counts, metadata = align_samples(
        counts,
        metadata,
        args.sample_id_column
    )

    # --------------------------------------------------------
    # 4. Confirm that counts are valid raw integer counts.
    # --------------------------------------------------------

    counts = validate_counts(counts)

    # Keep the original validated table for before/after QC statistics.
    original_counts = counts.copy()

    # --------------------------------------------------------
    # 5. Calculate sample-level QC.
    # --------------------------------------------------------

    sample_summary = make_sample_summary(
        counts,
        metadata,
        args.sample_id_column,
        args.group_column
    )

    # A zero-depth sample has no usable sequencing information.
    if (sample_summary["QC_Flag"] != "PASS").any():

        failed = sample_summary.loc[
            sample_summary["QC_Flag"] != "PASS",
            "SampleID"
        ].tolist()

        raise ValueError(
            "Zero-library sample(s) detected: "
            + ", ".join(failed)
        )

    # --------------------------------------------------------
    # 6. Calculate taxon-level QC.
    # --------------------------------------------------------

    taxon_summary = make_taxon_summary(
        counts,
        args.min_prevalence,
        args.min_total_count
    )

    # --------------------------------------------------------
    # 7. Apply the shared taxon filter.
    # --------------------------------------------------------

    counts = filter_taxa(
        counts,
        taxon_summary
    )

    # --------------------------------------------------------
    # 8. Build standardized output filenames.
    # --------------------------------------------------------

    prefix = f"qc_{dataset_label}"

    counts_file = (
        args.output_dir
        / f"{prefix}_counts.csv"
    )

    metadata_file = (
        args.output_dir
        / f"{prefix}_metadata.csv"
    )

    sample_file = (
        args.output_dir
        / f"{prefix}_sample_summary.csv"
    )

    taxon_file = (
        args.output_dir
        / f"{prefix}_taxon_summary.csv"
    )

    summary_file = (
        args.output_dir
        / f"{prefix}_summary.txt"
    )

    # --------------------------------------------------------
    # 9. Save the QC'd data used by all four methods.
    # --------------------------------------------------------

    counts.to_csv(
        counts_file
    )

    metadata.to_csv(
        metadata_file,
        index=False
    )

    # --------------------------------------------------------
    # 10. Save QC documentation.
    # --------------------------------------------------------

    sample_summary.to_csv(
        sample_file,
        index=False
    )

    taxon_summary.to_csv(
        taxon_file,
        index=False
    )

    write_summary(
        summary_file,
        dataset_label,
        original_counts,
        counts,
        metadata,
        sample_summary,
        taxon_summary,
        args.group_column,
        args.min_prevalence,
        args.min_total_count
    )

    # --------------------------------------------------------
    # 11. Show the pipeline handoff in the terminal.
    # --------------------------------------------------------

    print("QUALITY CONTROL COMPLETED")
    print("=========================")

    print(
        f"Samples: {counts.shape[1]}"
    )

    print(
        f"Taxa before QC: {original_counts.shape[0]}"
    )

    print(
        f"Taxa after QC: {counts.shape[0]}"
    )

    print(
        f"Taxa removed: "
        f"{original_counts.shape[0] - counts.shape[0]}"
    )

    print()
    print("Inputs for all four methods:")
    print(f"  {counts_file}")
    print(f"  {metadata_file}")

    print()
    print("QC reports:")
    print(f"  {sample_file}")
    print(f"  {taxon_file}")
    print(f"  {summary_file}")

    print()


# ============================================================
# RUN SCRIPT
# ============================================================

if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print()
        print("QUALITY CONTROL FAILED")
        print("======================")
        print(error)
        print()

        # A non-zero exit code lets a future pipeline know
        # that it should stop instead of running the DA methods.
        sys.exit(1)
