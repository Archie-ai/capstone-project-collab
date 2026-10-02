#!/usr/bin/env python3

# ============================================================
# RESULTS QC
# ============================================================
#
# This script checks that the standardized method results and
# downstream analysis outputs are internally consistent.
#
# No differential-abundance testing is performed here.
#
# This script was developed with assistance from OpenAI's ChatGPT
#

import argparse
from itertools import combinations
from pathlib import Path
import sys

import numpy as np
import pandas as pd


# ============================================================
# COMMAND-LINE INPUTS
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="QC the standardized and analyzed DA results."
    )

    parser.add_argument(
        "--standard-results",
        required=True,
        type=Path,
        help="Combined standardized results table."
    )

    parser.add_argument(
        "--analysis-dir",
        required=True,
        type=Path,
        help="Stage 5 analysis output folder."
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Folder for Stage 6 QC outputs."
    )

    parser.add_argument(
        "--dataset-label",
        required=True,
        help="Dataset label, such as null or real_test."
    )

    parser.add_argument(
        "--ground-truth",
        type=Path,
        default=None,
        help="Optional simulation ground-truth table."
    )

    return parser.parse_args()


# ============================================================
# BOOLEAN CLEANING
# ============================================================

def clean_boolean(series):
    """Convert common TRUE/FALSE and 1/0 values to booleans."""

    mapping = {
        True: True,
        False: False,
        1: True,
        0: False,
        "1": True,
        "0": False,
        "true": True,
        "false": False,
        "TRUE": True,
        "FALSE": False,
        "True": True,
        "False": False
    }

    cleaned = series.map(mapping)

    if cleaned.isna().any():
        raise ValueError(
            "A boolean column contains values that could not be interpreted."
        )

    return cleaned.astype(bool)


# ============================================================
# RECORD QC CHECK
# ============================================================

def add_check(checks, name, passed, details):
    """Add one PASS/FAIL result to the QC report."""

    checks.append({
        "Check": name,
        "Status": "PASS" if passed else "FAIL",
        "Details": details
    })


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    dataset = args.dataset_label

    if not args.standard_results.exists():
        raise FileNotFoundError(
            f"Standard results not found: {args.standard_results}"
        )

    if not args.analysis_dir.exists():
        raise FileNotFoundError(
            f"Analysis folder not found: {args.analysis_dir}"
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    checks = []

    print()
    print(f"Running results QC for dataset: {dataset}")
    print()


    # --------------------------------------------------------
    # 1. Load standardized results.
    # --------------------------------------------------------

    results = pd.read_csv(
        args.standard_results
    )

    required_columns = {
        "Dataset",
        "Method",
        "Taxon",
        "EffectSize",
        "EffectMetric",
        "PValue",
        "AdjustedPValue",
        "Significant",
        "Direction"
    }

    missing_columns = (
        required_columns - set(results.columns)
    )

    add_check(
        checks,
        "Required standardized columns",
        len(missing_columns) == 0,
        (
            "All required columns present"
            if not missing_columns
            else "Missing: " + ", ".join(sorted(missing_columns))
        )
    )

    if missing_columns:
        raise ValueError(
            "Standardized results are missing required columns."
        )

    results["Significant"] = clean_boolean(
        results["Significant"]
    )


    # --------------------------------------------------------
    # 2. Confirm all four methods are present.
    # --------------------------------------------------------

    expected_methods = {
        "CLR_Wilcoxon",
        "ALDEx2",
        "ANCOM_BC",
        "DESeq2"
    }

    observed_methods = set(
        results["Method"].unique()
    )

    add_check(
        checks,
        "All four methods present",
        observed_methods == expected_methods,
        (
            "Methods: "
            + ", ".join(sorted(observed_methods))
        )
    )


    # --------------------------------------------------------
    # 3. Check duplicate Method/Taxon rows.
    # --------------------------------------------------------

    duplicate_count = int(
        results.duplicated(
            subset=["Method", "Taxon"]
        ).sum()
    )

    add_check(
        checks,
        "No duplicate Method/Taxon rows",
        duplicate_count == 0,
        f"Duplicate rows: {duplicate_count}"
    )


    # --------------------------------------------------------
    # 4. Check that each method analyzed the same taxa.
    # --------------------------------------------------------

    taxa_by_method = {
        method: set(group["Taxon"])
        for method, group in results.groupby("Method")
    }

    taxa_sets = list(
        taxa_by_method.values()
    )

    same_taxa = all(
        taxa_set == taxa_sets[0]
        for taxa_set in taxa_sets[1:]
    )

    counts_text = ", ".join(
        f"{method}={len(taxa)}"
        for method, taxa in sorted(taxa_by_method.items())
    )

    add_check(
        checks,
        "Same taxa across methods",
        same_taxa,
        counts_text
    )


    # --------------------------------------------------------
    # 5. Check p-value ranges.
    # --------------------------------------------------------

    p_valid = (
        results["PValue"]
        .dropna()
        .between(0, 1)
        .all()
    )

    adjusted_valid = (
        results["AdjustedPValue"]
        .dropna()
        .between(0, 1)
        .all()
    )

    add_check(
        checks,
        "Raw p-values within 0 to 1",
        bool(p_valid),
        "Checked all non-missing raw p-values"
    )

    add_check(
        checks,
        "Adjusted p-values within 0 to 1",
        bool(adjusted_valid),
        "Checked all non-missing adjusted p-values"
    )


    # --------------------------------------------------------
    # 6. Confirm significance matches adjusted p-value < 0.05.
    # --------------------------------------------------------

    expected_significant = (
        results["AdjustedPValue"].notna()
        & (results["AdjustedPValue"] < 0.05)
    )

    significance_matches = (
        results["Significant"]
        == expected_significant
    )

    mismatch_count = int(
        (~significance_matches).sum()
    )

    add_check(
        checks,
        "Significance matches adjusted p-value threshold",
        mismatch_count == 0,
        f"Mismatched rows: {mismatch_count}"
    )


    # --------------------------------------------------------
    # 7. Confirm direction matches effect-size sign.
    # --------------------------------------------------------

    direction_mismatches = 0

    for _, row in results.iterrows():

        effect = row["EffectSize"]
        direction = row["Direction"]

        if pd.isna(effect):
            continue

        if effect > 0:
            expected_direction = "Higher_in_Case"

        elif effect < 0:
            expected_direction = "Higher_in_Control"

        else:
            expected_direction = "No_change"

        if direction != expected_direction:
            direction_mismatches += 1

    add_check(
        checks,
        "Effect direction matches effect-size sign",
        direction_mismatches == 0,
        f"Mismatched rows: {direction_mismatches}"
    )


    # --------------------------------------------------------
    # 8. Locate Stage 5 analysis files.
    # --------------------------------------------------------

    significant_file = (
        args.analysis_dir
        / f"significant_calls_{dataset}.csv"
    )

    consensus_file = (
        args.analysis_dir
        / f"taxon_consensus_{dataset}.csv"
    )

    overlap_file = (
        args.analysis_dir
        / f"pairwise_overlap_{dataset}.csv"
    )

    analysis_summary_file = (
        args.analysis_dir
        / f"analysis_{dataset}_summary.txt"
    )

    analysis_files = [
        significant_file,
        consensus_file,
        overlap_file,
        analysis_summary_file
    ]

    missing_analysis_files = [
        str(file)
        for file in analysis_files
        if not file.exists()
    ]

    add_check(
        checks,
        "Required Stage 5 analysis files present",
        len(missing_analysis_files) == 0,
        (
            "All required files present"
            if not missing_analysis_files
            else "Missing: " + "; ".join(missing_analysis_files)
        )
    )

    if missing_analysis_files:
        raise FileNotFoundError(
            "One or more required Stage 5 files are missing."
        )


    # --------------------------------------------------------
    # 9. Check significant_calls file against standardized data.
    # --------------------------------------------------------

    significant_calls = pd.read_csv(
        significant_file
    )

    expected_calls = results[
        results["Significant"]
    ][["Method", "Taxon"]].copy()

    observed_calls = significant_calls[
        ["Method", "Taxon"]
    ].copy()

    expected_call_set = set(
        map(tuple, expected_calls.to_numpy())
    )

    observed_call_set = set(
        map(tuple, observed_calls.to_numpy())
    )

    calls_match = (
        expected_call_set == observed_call_set
    )

    add_check(
        checks,
        "Stage 5 significant calls match standardized results",
        calls_match,
        (
            f"Expected calls: {len(expected_call_set)}; "
            f"Observed calls: {len(observed_call_set)}"
        )
    )


    # --------------------------------------------------------
    # 10. Check taxon consensus counts.
    # --------------------------------------------------------

    consensus = pd.read_csv(
        consensus_file
    )

    expected_consensus = (
        results
        .groupby("Taxon")["Significant"]
        .sum()
        .astype(int)
    )

    consensus_lookup = (
        consensus
        .set_index("Taxon")["SignificantMethodCount"]
    )

    common_taxa = (
        expected_consensus.index
        .intersection(consensus_lookup.index)
    )

    consensus_count_matches = (
        expected_consensus.loc[common_taxa]
        == consensus_lookup.loc[common_taxa]
    )

    consensus_mismatches = int(
        (~consensus_count_matches).sum()
    )

    same_consensus_taxa = (
        set(expected_consensus.index)
        == set(consensus_lookup.index)
    )

    add_check(
        checks,
        "Taxon consensus counts match method calls",
        (
            same_consensus_taxa
            and consensus_mismatches == 0
        ),
        (
            f"Consensus count mismatches: {consensus_mismatches}"
        )
    )


    # --------------------------------------------------------
    # 11. Recheck pairwise significant overlap.
    # --------------------------------------------------------

    overlap = pd.read_csv(
        overlap_file
    )

    overlap_errors = 0

    methods = sorted(
        results["Method"].unique()
    )

    for method_a, method_b in combinations(methods, 2):

        set_a = set(
            results.loc[
                (results["Method"] == method_a)
                & results["Significant"],
                "Taxon"
            ]
        )

        set_b = set(
            results.loc[
                (results["Method"] == method_b)
                & results["Significant"],
                "Taxon"
            ]
        )

        expected_shared = len(
            set_a & set_b
        )

        row = overlap[
            (
                (overlap["MethodA"] == method_a)
                & (overlap["MethodB"] == method_b)
            )
        ]

        if len(row) != 1:
            overlap_errors += 1
            continue

        observed_shared = int(
            row.iloc[0]["SharedSignificant"]
        )

        if observed_shared != expected_shared:
            overlap_errors += 1

    add_check(
        checks,
        "Pairwise overlap matches underlying calls",
        overlap_errors == 0,
        f"Pairwise overlap errors: {overlap_errors}"
    )


    # --------------------------------------------------------
    # 12. Check ground-truth benchmark when supplied.
    # --------------------------------------------------------

    if args.ground_truth is not None:

        if not args.ground_truth.exists():
            raise FileNotFoundError(
                f"Ground truth not found: {args.ground_truth}"
            )

        performance_file = (
            args.analysis_dir
            / f"method_performance_{dataset}.csv"
        )

        if not performance_file.exists():
            raise FileNotFoundError(
                f"Performance table not found: {performance_file}"
            )

        truth = pd.read_csv(
            args.ground_truth
        )

        truth["True_DA"] = clean_boolean(
            truth["True_DA"]
        )

        performance = pd.read_csv(
            performance_file
        )

        merged = results.merge(
            truth[["Taxon", "True_DA"]],
            on="Taxon",
            how="left",
            validate="many_to_one"
        )

        performance_errors = 0

        for method, method_data in merged.groupby("Method"):

            predicted = (
                method_data["Significant"]
            )

            true_da = (
                method_data["True_DA"]
            )

            tp = int(
                (predicted & true_da).sum()
            )

            fp = int(
                (predicted & ~true_da).sum()
            )

            tn = int(
                (~predicted & ~true_da).sum()
            )

            fn = int(
                (~predicted & true_da).sum()
            )

            row = performance[
                performance["Method"] == method
            ]

            if len(row) != 1:
                performance_errors += 1
                continue

            row = row.iloc[0]

            if (
                int(row["TruePositive"]) != tp
                or int(row["FalsePositive"]) != fp
                or int(row["TrueNegative"]) != tn
                or int(row["FalseNegative"]) != fn
            ):
                performance_errors += 1

        add_check(
            checks,
            "Ground-truth performance counts are reproducible",
            performance_errors == 0,
            f"Method performance errors: {performance_errors}"
        )


    # --------------------------------------------------------
    # 13. Save detailed QC table.
    # --------------------------------------------------------

    checks_df = pd.DataFrame(
        checks
    )

    checks_file = (
        args.output_dir
        / f"results_qc_{dataset}_checks.csv"
    )

    checks_df.to_csv(
        checks_file,
        index=False
    )


    # --------------------------------------------------------
    # 14. Decide overall PASS / FAIL.
    # --------------------------------------------------------

    failed_checks = checks_df[
        checks_df["Status"] == "FAIL"
    ]

    overall_status = (
        "PASS"
        if failed_checks.empty
        else "FAIL"
    )


    # --------------------------------------------------------
    # 15. Write readable QC summary.
    # --------------------------------------------------------

    summary_file = (
        args.output_dir
        / f"results_qc_{dataset}_summary.txt"
    )

    lines = [
        "RESULTS QUALITY-CONTROL SUMMARY",
        "===============================",
        "",
        f"Dataset: {dataset}",
        f"Overall status: {overall_status}",
        "",
        f"Checks performed: {len(checks_df)}",
        f"Checks passed: {(checks_df['Status'] == 'PASS').sum()}",
        f"Checks failed: {(checks_df['Status'] == 'FAIL').sum()}",
        "",
        "CHECK RESULTS",
        "-------------"
    ]

    for _, row in checks_df.iterrows():

        lines.append(
            f"{row['Status']} - {row['Check']}: "
            f"{row['Details']}"
        )

    summary_file.write_text(
        "\n".join(lines) + "\n"
    )


    # --------------------------------------------------------
    # 16. Terminal output.
    # --------------------------------------------------------

    print("RESULTS QC COMPLETED")
    print("====================")

    print(
        f"Overall status: {overall_status}"
    )

    print(
        f"Checks passed: "
        f"{(checks_df['Status'] == 'PASS').sum()}"
    )

    print(
        f"Checks failed: "
        f"{(checks_df['Status'] == 'FAIL').sum()}"
    )

    print()
    print("Detailed checks:")
    print(f"  {checks_file}")

    print()
    print("Summary:")
    print(f"  {summary_file}")

    print()


    # A failed QC stage should stop the future automated pipeline.
    if overall_status == "FAIL":
        sys.exit(1)


# ============================================================
# RUN SCRIPT
# ============================================================

if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print()
        print("RESULTS QC FAILED")
        print("=================")
        print(error)
        print()

        sys.exit(1)
