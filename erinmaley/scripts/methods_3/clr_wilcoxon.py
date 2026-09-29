#!/usr/bin/env python3

# ============================================================
# CLR / WILCOXON
# ============================================================
#
# This script runs CLR transformation followed by a Wilcoxon
# rank-sum test on microbiome count data that passed shared QC.
#
# This script was developed with assistance from OpenAI's ChatGPT
#

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


# ============================================================
# COMMAND-LINE INPUTS
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="Run CLR/Wilcoxon differential abundance analysis."
    )

    parser.add_argument(
        "--counts",
        required=True,
        type=Path,
        help="QC count table."
    )

    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="QC metadata table."
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Folder for CLR/Wilcoxon outputs."
    )

    parser.add_argument(
        "--dataset-label",
        required=True,
        help="Dataset label, e.g. null or real_test."
    )

    return parser.parse_args()


# ============================================================
# BENJAMINI-HOCHBERG CORRECTION
# ============================================================

def bh_adjust(p_values):
    """
    Adjust p-values using the Benjamini-Hochberg FDR procedure.

    This keeps the significance criterion consistent with the
    other differential-abundance methods.
    """

    p_values = np.asarray(p_values, dtype=float)

    adjusted = np.full(
        len(p_values),
        np.nan
    )

    valid = ~np.isnan(p_values)

    p = p_values[valid]

    if len(p) == 0:
        return adjusted

    order = np.argsort(p)

    ranked_p = p[order]

    n = len(ranked_p)

    adjusted_ranked = (
        ranked_p * n / np.arange(1, n + 1)
    )

    # BH adjusted p-values must not decrease as rank increases.
    adjusted_ranked = np.minimum.accumulate(
        adjusted_ranked[::-1]
    )[::-1]

    adjusted_ranked = np.minimum(
        adjusted_ranked,
        1.0
    )

    restored = np.empty(n)

    restored[order] = adjusted_ranked

    adjusted[valid] = restored

    return adjusted


# ============================================================
# CLR TRANSFORMATION
# ============================================================

def clr_transform(counts, pseudocount=1.0):
    """
    Apply centered log-ratio transformation sample by sample.

    A pseudocount of 1 is added because log(0) is undefined.

    For each sample:
        CLR abundance =
        log(taxon abundance) - mean(log(all taxa in sample))
    """

    adjusted_counts = counts + pseudocount

    log_counts = np.log(
        adjusted_counts
    )

    clr = log_counts.subtract(
        log_counts.mean(axis=0),
        axis=1
    )

    return clr


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():

    args = parse_args()

    if not args.counts.exists():
        raise FileNotFoundError(
            f"Count file not found: {args.counts}"
        )

    if not args.metadata.exists():
        raise FileNotFoundError(
            f"Metadata file not found: {args.metadata}"
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # 1. Load QC-approved input.
    # --------------------------------------------------------

    counts = pd.read_csv(
        args.counts,
        index_col=0
    )

    metadata = pd.read_csv(
        args.metadata
    )

    print()
    print(
        f"Running CLR/Wilcoxon for dataset: "
        f"{args.dataset_label}"
    )
    print()

    print(f"Input taxa: {counts.shape[0]}")
    print(f"Input samples: {counts.shape[1]}")


    # --------------------------------------------------------
    # 2. Check metadata and sample matching.
    # --------------------------------------------------------

    required_columns = {
        "SampleID",
        "Group"
    }

    missing = (
        required_columns
        - set(metadata.columns)
    )

    if missing:
        raise ValueError(
            "Metadata is missing required column(s): "
            + ", ".join(sorted(missing))
        )

    metadata["SampleID"] = (
        metadata["SampleID"]
        .astype(str)
    )

    counts.columns = (
        counts.columns
        .astype(str)
    )

    if set(counts.columns) != set(metadata["SampleID"]):
        raise ValueError(
            "Sample IDs do not match between counts and metadata."
        )

    metadata = (
        metadata
        .set_index("SampleID")
        .loc[counts.columns]
    )

    control_samples = metadata.index[
        metadata["Group"] == "Control"
    ].tolist()

    case_samples = metadata.index[
        metadata["Group"] == "Case"
    ].tolist()

    if not control_samples or not case_samples:
        raise ValueError(
            "Both Control and Case samples are required."
        )

    print(
        f"Groups: Control {len(control_samples)}, "
        f"Case {len(case_samples)}"
    )


    # --------------------------------------------------------
    # 3. CLR transform the counts.
    # --------------------------------------------------------

    # Zero handling is method-specific, so it happens here
    # rather than during shared QC.
    pseudocount = 1.0

    clr = clr_transform(
        counts,
        pseudocount=pseudocount
    )


    # --------------------------------------------------------
    # 4. Run Wilcoxon rank-sum test for every taxon.
    # --------------------------------------------------------

    rows = []

    for taxon in clr.index:

        control_values = clr.loc[
            taxon,
            control_samples
        ].to_numpy()

        case_values = clr.loc[
            taxon,
            case_samples
        ].to_numpy()

        # Mann-Whitney U is the Wilcoxon rank-sum test for
        # two independent groups.
        test = mannwhitneyu(
            case_values,
            control_values,
            alternative="two-sided",
            method="auto"
        )

        control_median = np.median(
            control_values
        )

        case_median = np.median(
            case_values
        )

        # Positive difference means higher CLR abundance in Case.
        clr_difference = (
            case_median
            - control_median
        )

        rows.append({
            "Taxon": taxon,
            "CaseMedianCLR": case_median,
            "ControlMedianCLR": control_median,
            "CLRDifference": clr_difference,
            "Statistic": test.statistic,
            "PValue": test.pvalue
        })


    results = pd.DataFrame(rows)


    # --------------------------------------------------------
    # 5. Correct for multiple testing.
    # --------------------------------------------------------

    results["AdjustedPValue"] = bh_adjust(
        results["PValue"]
    )

    results["Significant"] = (
        results["AdjustedPValue"] < 0.05
    )


    # --------------------------------------------------------
    # 6. Add direction.
    # --------------------------------------------------------

    results["Direction"] = np.where(
        results["CLRDifference"] > 0,
        "Higher_in_Case",
        np.where(
            results["CLRDifference"] < 0,
            "Higher_in_Control",
            "No_change"
        )
    )


    # --------------------------------------------------------
    # 7. Save result table.
    # --------------------------------------------------------

    output_file = (
        args.output_dir
        / f"clr_wilcoxon_{args.dataset_label}_results.csv"
    )

    results.to_csv(
        output_file,
        index=False
    )


    # --------------------------------------------------------
    # 8. Save run summary.
    # --------------------------------------------------------

    significant_count = int(
        results["Significant"].sum()
    )

    summary_file = (
        args.output_dir
        / f"clr_wilcoxon_{args.dataset_label}_summary.txt"
    )

    summary = f"""
CLR / WILCOXON METHOD SUMMARY
=============================

Dataset:
{args.dataset_label}

Samples:
{counts.shape[1]}

Taxa analyzed:
{counts.shape[0]}

Comparison:
Case vs Control

Transformation:
Centered log-ratio (CLR)

Zero handling:
Pseudocount = {pseudocount}

Statistical test:
Two-sided Wilcoxon rank-sum / Mann-Whitney U

P-value adjustment:
Benjamini-Hochberg (BH)

Significance threshold:
Adjusted p-value < 0.05

Significant taxa:
{significant_count}

Direction:
Positive CLR difference = higher in Case
Negative CLR difference = higher in Control

QC filtering was completed before CLR/Wilcoxon.
"""

    summary_file.write_text(
        summary.strip() + "\n"
    )


    # --------------------------------------------------------
    # 9. Terminal output.
    # --------------------------------------------------------

    print()
    print("CLR/WILCOXON COMPLETED")
    print("=======================")

    print(
        f"Taxa analyzed: {len(results)}"
    )

    print(
        f"Significant taxa: {significant_count}"
    )

    print()
    print("Output:")
    print(f"  {output_file}")

    print()
    print("Summary:")
    print(f"  {summary_file}")

    print()


if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print()
        print("CLR/WILCOXON FAILED")
        print("====================")
        print(error)
        print()

        sys.exit(1)
