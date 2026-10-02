#!/usr/bin/env python3

# ============================================================
# ANALYSIS
# ============================================================
#
# This script compares the standardized results from all four
# differential-abundance methods.
#
# If simulation ground truth is supplied, it also evaluates
# true/false positive performance.
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
        description="Analyze standardized differential-abundance results."
    )

    parser.add_argument(
        "--standard-results",
        required=True,
        type=Path,
        help="Combined standardized results table."
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Folder where analysis outputs will be saved."
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
    """Convert common TRUE/FALSE and 1/0 values into Python booleans."""

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
# SAFE DIVISION
# ============================================================

def safe_divide(numerator, denominator):
    """Return NA when a performance metric has no valid denominator."""

    if denominator == 0:
        return np.nan

    return numerator / denominator


# ============================================================
# METHOD PERFORMANCE
# ============================================================

def calculate_method_performance(results):
    """
    Compare each method's significance calls against known truth.

    TP = significant and truly DA
    FP = significant but not truly DA
    TN = not significant and not truly DA
    FN = not significant but truly DA
    """

    rows = []

    for method, method_data in results.groupby("Method"):

        predicted = method_data["Significant"]
        truth = method_data["True_DA"]

        tp = int((predicted & truth).sum())
        fp = int((predicted & ~truth).sum())
        tn = int((~predicted & ~truth).sum())
        fn = int((~predicted & truth).sum())

        rows.append({
            "Method": method,
            "TaxaTested": len(method_data),
            "SignificantTaxa": int(predicted.sum()),
            "TruePositive": tp,
            "FalsePositive": fp,
            "TrueNegative": tn,
            "FalseNegative": fn,

            "Sensitivity": safe_divide(
                tp,
                tp + fn
            ),

            "Specificity": safe_divide(
                tn,
                tn + fp
            ),

            "Precision": safe_divide(
                tp,
                tp + fp
            ),

            "FalseDiscoveryRate": safe_divide(
                fp,
                tp + fp
            ),

            "FalsePositiveRate": safe_divide(
                fp,
                fp + tn
            ),

            "Accuracy": safe_divide(
                tp + tn,
                tp + fp + tn + fn
            )
        })

    return pd.DataFrame(rows)


# ============================================================
# PAIRWISE METHOD OVERLAP
# ============================================================

def calculate_pairwise_overlap(results):
    """Measure overlap between significant-taxon sets for each method pair."""

    methods = sorted(
        results["Method"].unique()
    )

    rows = []

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

        intersection = set_a & set_b
        union = set_a | set_b

        jaccard = (
            len(intersection) / len(union)
            if len(union) > 0
            else np.nan
        )

        rows.append({
            "MethodA": method_a,
            "MethodB": method_b,
            "SignificantA": len(set_a),
            "SignificantB": len(set_b),
            "SharedSignificant": len(intersection),
            "UnionSignificant": len(union),
            "JaccardIndex": jaccard,
            "SharedTaxa": ";".join(sorted(intersection))
        })

    return pd.DataFrame(rows)


# ============================================================
# TAXON-LEVEL CONSENSUS
# ============================================================

def calculate_taxon_consensus(results):
    """
    Summarize how many methods called each taxon significant.

    This table becomes useful later for the consensus/report stage.
    """

    rows = []

    method_count = results["Method"].nunique()

    for taxon, taxon_data in results.groupby("Taxon"):

        significant_data = taxon_data[
            taxon_data["Significant"]
        ]

        methods_significant = sorted(
            significant_data["Method"].tolist()
        )

        directions = sorted(
            significant_data["Direction"]
            .dropna()
            .unique()
            .tolist()
        )

        significant_count = len(
            methods_significant
        )

        if significant_count == 0:
            direction_agreement = "No_significant_calls"

        elif significant_count == 1:
            direction_agreement = "Not_applicable"

        elif len(directions) == 1:
            direction_agreement = "Agree"

        else:
            direction_agreement = "Disagree"

        row = {
            "Taxon": taxon,
            "MethodsSignificant": ";".join(methods_significant),
            "SignificantMethodCount": significant_count,
            "TotalMethods": method_count,
            "ConsensusFraction": significant_count / method_count,
            "SignificantDirections": ";".join(directions),
            "DirectionAgreement": direction_agreement
        }

        # Include ground truth when available.
        if "True_DA" in taxon_data.columns:

            row["True_DA"] = bool(
                taxon_data["True_DA"].iloc[0]
            )

            row["True_Direction"] = (
                taxon_data["True_Direction"].iloc[0]
            )

        rows.append(row)

    consensus = pd.DataFrame(rows)

    return consensus.sort_values(
        by=[
            "SignificantMethodCount",
            "Taxon"
        ],
        ascending=[
            False,
            True
        ]
    )


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    if not args.standard_results.exists():
        raise FileNotFoundError(
            f"Standardized result file not found: "
            f"{args.standard_results}"
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    dataset = args.dataset_label

    print()
    print(f"Analyzing dataset: {dataset}")
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

    missing = required_columns - set(
        results.columns
    )

    if missing:
        raise ValueError(
            "Standardized results are missing column(s): "
            + ", ".join(sorted(missing))
        )

    results["Significant"] = clean_boolean(
        results["Significant"]
    )

    print(
        f"Methods: {results['Method'].nunique()}"
    )

    print(
        f"Unique taxa: {results['Taxon'].nunique()}"
    )

    print(
        f"Total method-taxon results: {len(results)}"
    )


    # --------------------------------------------------------
    # 2. Check for duplicate method/taxon rows.
    # --------------------------------------------------------

    duplicates = results.duplicated(
        subset=[
            "Method",
            "Taxon"
        ]
    )

    if duplicates.any():
        raise ValueError(
            "Duplicate Method/Taxon rows found in standardized results."
        )


    # --------------------------------------------------------
    # 3. Add simulation ground truth when supplied.
    # --------------------------------------------------------

    truth_available = (
        args.ground_truth is not None
    )

    if truth_available:

        if not args.ground_truth.exists():
            raise FileNotFoundError(
                f"Ground-truth file not found: "
                f"{args.ground_truth}"
            )

        truth = pd.read_csv(
            args.ground_truth
        )

        required_truth = {
            "Taxon",
            "True_DA",
            "True_Direction"
        }

        missing_truth = (
            required_truth - set(truth.columns)
        )

        if missing_truth:
            raise ValueError(
                "Ground truth is missing column(s): "
                + ", ".join(sorted(missing_truth))
            )

        truth["True_DA"] = clean_boolean(
            truth["True_DA"]
        )

        results = results.merge(
            truth,
            on="Taxon",
            how="left",
            validate="many_to_one"
        )

        if results["True_DA"].isna().any():
            raise ValueError(
                "Some analyzed taxa are missing from ground truth."
            )

        print(
            f"Ground-truth DA taxa: "
            f"{int(truth['True_DA'].sum())}"
        )


    # --------------------------------------------------------
    # 4. Save significant calls from all methods.
    # --------------------------------------------------------

    significant_calls = results[
        results["Significant"]
    ].copy()

    significant_file = (
        args.output_dir
        / f"significant_calls_{dataset}.csv"
    )

    significant_calls.to_csv(
        significant_file,
        index=False
    )


    # --------------------------------------------------------
    # 5. Calculate pairwise overlap.
    # --------------------------------------------------------

    overlap = calculate_pairwise_overlap(
        results
    )

    overlap_file = (
        args.output_dir
        / f"pairwise_overlap_{dataset}.csv"
    )

    overlap.to_csv(
        overlap_file,
        index=False
    )


    # --------------------------------------------------------
    # 6. Calculate taxon-level consensus.
    # --------------------------------------------------------

    consensus = calculate_taxon_consensus(
        results
    )

    consensus_file = (
        args.output_dir
        / f"taxon_consensus_{dataset}.csv"
    )

    consensus.to_csv(
        consensus_file,
        index=False
    )


    # --------------------------------------------------------
    # 7. Evaluate benchmark performance when truth exists.
    # --------------------------------------------------------

    performance = None

    if truth_available:

        performance = calculate_method_performance(
            results
        )

        performance_file = (
            args.output_dir
            / f"method_performance_{dataset}.csv"
        )

        performance.to_csv(
            performance_file,
            index=False
        )


    # --------------------------------------------------------
    # 8. Write readable analysis summary.
    # --------------------------------------------------------

    summary_file = (
        args.output_dir
        / f"analysis_{dataset}_summary.txt"
    )

    method_calls = (
        results
        .groupby("Method")["Significant"]
        .sum()
        .sort_index()
    )

    lines = [
        "DIFFERENTIAL-ABUNDANCE ANALYSIS SUMMARY",
        "=======================================",
        "",
        f"Dataset: {dataset}",
        f"Methods: {results['Method'].nunique()}",
        f"Unique taxa: {results['Taxon'].nunique()}",
        "",
        "SIGNIFICANT TAXA BY METHOD",
        "--------------------------"
    ]

    for method, count in method_calls.items():

        lines.append(
            f"{method}: {int(count)}"
        )

    lines.extend([
        "",
        "TAXON CONSENSUS",
        "---------------",
        (
            "Taxa called significant by at least one method: "
            f"{int((consensus['SignificantMethodCount'] > 0).sum())}"
        ),
        (
            "Taxa called significant by all methods: "
            f"{int((consensus['SignificantMethodCount'] == results['Method'].nunique()).sum())}"
        )
    ])


    # Ground-truth performance is only meaningful for simulated
    # datasets where the true DA status is known.
    if truth_available:

        lines.extend([
            "",
            "GROUND-TRUTH BENCHMARK",
            "----------------------",
            f"True DA taxa: {int(results.drop_duplicates('Taxon')['True_DA'].sum())}",
            ""
        ])

        for _, row in performance.iterrows():

            lines.extend([
                row["Method"],
                f"  True positives: {row['TruePositive']}",
                f"  False positives: {row['FalsePositive']}",
                f"  True negatives: {row['TrueNegative']}",
                f"  False negatives: {row['FalseNegative']}",
                f"  Specificity: {row['Specificity']:.4f}"
                if not pd.isna(row["Specificity"])
                else "  Specificity: NA",
                f"  False positive rate: {row['FalsePositiveRate']:.4f}"
                if not pd.isna(row["FalsePositiveRate"])
                else "  False positive rate: NA",
                ""
            ])


    summary_file.write_text(
        "\n".join(lines) + "\n"
    )


    # --------------------------------------------------------
    # 9. Terminal summary.
    # --------------------------------------------------------

    print()
    print("ANALYSIS COMPLETED")
    print("==================")

    print()
    print("Significant taxa by method:")

    for method, count in method_calls.items():
        print(
            f"  {method}: {int(count)}"
        )

    print()

    if truth_available:
        print(
            "Ground-truth performance table:"
        )
        print(
            f"  {performance_file}"
        )

    print(
        "Taxon consensus table:"
    )
    print(
        f"  {consensus_file}"
    )

    print(
        "Pairwise overlap table:"
    )
    print(
        f"  {overlap_file}"
    )

    print(
        "Analysis summary:"
    )
    print(
        f"  {summary_file}"
    )

    print()


# ============================================================
# RUN SCRIPT
# ============================================================

if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print()
        print("ANALYSIS FAILED")
        print("===============")
        print(error)
        print()

        sys.exit(1)#!/usr/bin/env python3

# ============================================================
# ANALYSIS
# ============================================================
#
# This script compares the standardized results from all four
# differential-abundance methods.
#
# If simulation ground truth is supplied, it also evaluates
# true/false positive performance.
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
        description="Analyze standardized differential-abundance results."
    )

    parser.add_argument(
        "--standard-results",
        required=True,
        type=Path,
        help="Combined standardized results table."
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Folder where analysis outputs will be saved."
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
    """Convert common TRUE/FALSE and 1/0 values into Python booleans."""

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
# SAFE DIVISION
# ============================================================

def safe_divide(numerator, denominator):
    """Return NA when a performance metric has no valid denominator."""

    if denominator == 0:
        return np.nan

    return numerator / denominator


# ============================================================
# METHOD PERFORMANCE
# ============================================================

def calculate_method_performance(results):
    """
    Compare each method's significance calls against known truth.

    TP = significant and truly DA
    FP = significant but not truly DA
    TN = not significant and not truly DA
    FN = not significant but truly DA
    """

    rows = []

    for method, method_data in results.groupby("Method"):

        predicted = method_data["Significant"]
        truth = method_data["True_DA"]

        tp = int((predicted & truth).sum())
        fp = int((predicted & ~truth).sum())
        tn = int((~predicted & ~truth).sum())
        fn = int((~predicted & truth).sum())

        rows.append({
            "Method": method,
            "TaxaTested": len(method_data),
            "SignificantTaxa": int(predicted.sum()),
            "TruePositive": tp,
            "FalsePositive": fp,
            "TrueNegative": tn,
            "FalseNegative": fn,

            "Sensitivity": safe_divide(
                tp,
                tp + fn
            ),

            "Specificity": safe_divide(
                tn,
                tn + fp
            ),

            "Precision": safe_divide(
                tp,
                tp + fp
            ),

            "FalseDiscoveryRate": safe_divide(
                fp,
                tp + fp
            ),

            "FalsePositiveRate": safe_divide(
                fp,
                fp + tn
            ),

            "Accuracy": safe_divide(
                tp + tn,
                tp + fp + tn + fn
            )
        })

    return pd.DataFrame(rows)


# ============================================================
# PAIRWISE METHOD OVERLAP
# ============================================================

def calculate_pairwise_overlap(results):
    """Measure overlap between significant-taxon sets for each method pair."""

    methods = sorted(
        results["Method"].unique()
    )

    rows = []

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

        intersection = set_a & set_b
        union = set_a | set_b

        jaccard = (
            len(intersection) / len(union)
            if len(union) > 0
            else np.nan
        )

        rows.append({
            "MethodA": method_a,
            "MethodB": method_b,
            "SignificantA": len(set_a),
            "SignificantB": len(set_b),
            "SharedSignificant": len(intersection),
            "UnionSignificant": len(union),
            "JaccardIndex": jaccard,
            "SharedTaxa": ";".join(sorted(intersection))
        })

    return pd.DataFrame(rows)


# ============================================================
# TAXON-LEVEL CONSENSUS
# ============================================================

def calculate_taxon_consensus(results):
    """
    Summarize how many methods called each taxon significant.

    This table becomes useful later for the consensus/report stage.
    """

    rows = []

    method_count = results["Method"].nunique()

    for taxon, taxon_data in results.groupby("Taxon"):

        significant_data = taxon_data[
            taxon_data["Significant"]
        ]

        methods_significant = sorted(
            significant_data["Method"].tolist()
        )

        directions = sorted(
            significant_data["Direction"]
            .dropna()
            .unique()
            .tolist()
        )

        significant_count = len(
            methods_significant
        )

        if significant_count == 0:
            direction_agreement = "No_significant_calls"

        elif significant_count == 1:
            direction_agreement = "Not_applicable"

        elif len(directions) == 1:
            direction_agreement = "Agree"

        else:
            direction_agreement = "Disagree"

        row = {
            "Taxon": taxon,
            "MethodsSignificant": ";".join(methods_significant),
            "SignificantMethodCount": significant_count,
            "TotalMethods": method_count,
            "ConsensusFraction": significant_count / method_count,
            "SignificantDirections": ";".join(directions),
            "DirectionAgreement": direction_agreement
        }

        # Include ground truth when available.
        if "True_DA" in taxon_data.columns:

            row["True_DA"] = bool(
                taxon_data["True_DA"].iloc[0]
            )

            row["True_Direction"] = (
                taxon_data["True_Direction"].iloc[0]
            )

        rows.append(row)

    consensus = pd.DataFrame(rows)

    return consensus.sort_values(
        by=[
            "SignificantMethodCount",
            "Taxon"
        ],
        ascending=[
            False,
            True
        ]
    )


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    if not args.standard_results.exists():
        raise FileNotFoundError(
            f"Standardized result file not found: "
            f"{args.standard_results}"
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    dataset = args.dataset_label

    print()
    print(f"Analyzing dataset: {dataset}")
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

    missing = required_columns - set(
        results.columns
    )

    if missing:
        raise ValueError(
            "Standardized results are missing column(s): "
            + ", ".join(sorted(missing))
        )

    results["Significant"] = clean_boolean(
        results["Significant"]
    )

    print(
        f"Methods: {results['Method'].nunique()}"
    )

    print(
        f"Unique taxa: {results['Taxon'].nunique()}"
    )

    print(
        f"Total method-taxon results: {len(results)}"
    )


    # --------------------------------------------------------
    # 2. Check for duplicate method/taxon rows.
    # --------------------------------------------------------

    duplicates = results.duplicated(
        subset=[
            "Method",
            "Taxon"
        ]
    )

    if duplicates.any():
        raise ValueError(
            "Duplicate Method/Taxon rows found in standardized results."
        )


    # --------------------------------------------------------
    # 3. Add simulation ground truth when supplied.
    # --------------------------------------------------------

    truth_available = (
        args.ground_truth is not None
    )

    if truth_available:

        if not args.ground_truth.exists():
            raise FileNotFoundError(
                f"Ground-truth file not found: "
                f"{args.ground_truth}"
            )

        truth = pd.read_csv(
            args.ground_truth
        )

        required_truth = {
            "Taxon",
            "True_DA",
            "True_Direction"
        }

        missing_truth = (
            required_truth - set(truth.columns)
        )

        if missing_truth:
            raise ValueError(
                "Ground truth is missing column(s): "
                + ", ".join(sorted(missing_truth))
            )

        truth["True_DA"] = clean_boolean(
            truth["True_DA"]
        )

        results = results.merge(
            truth,
            on="Taxon",
            how="left",
            validate="many_to_one"
        )

        if results["True_DA"].isna().any():
            raise ValueError(
                "Some analyzed taxa are missing from ground truth."
            )

        print(
            f"Ground-truth DA taxa: "
            f"{int(truth['True_DA'].sum())}"
        )


    # --------------------------------------------------------
    # 4. Save significant calls from all methods.
    # --------------------------------------------------------

    significant_calls = results[
        results["Significant"]
    ].copy()

    significant_file = (
        args.output_dir
        / f"significant_calls_{dataset}.csv"
    )

    significant_calls.to_csv(
        significant_file,
        index=False
    )


    # --------------------------------------------------------
    # 5. Calculate pairwise overlap.
    # --------------------------------------------------------

    overlap = calculate_pairwise_overlap(
        results
    )

    overlap_file = (
        args.output_dir
        / f"pairwise_overlap_{dataset}.csv"
    )

    overlap.to_csv(
        overlap_file,
        index=False
    )


    # --------------------------------------------------------
    # 6. Calculate taxon-level consensus.
    # --------------------------------------------------------

    consensus = calculate_taxon_consensus(
        results
    )

    consensus_file = (
        args.output_dir
        / f"taxon_consensus_{dataset}.csv"
    )

    consensus.to_csv(
        consensus_file,
        index=False
    )


    # --------------------------------------------------------
    # 7. Evaluate benchmark performance when truth exists.
    # --------------------------------------------------------

    performance = None

    if truth_available:

        performance = calculate_method_performance(
            results
        )

        performance_file = (
            args.output_dir
            / f"method_performance_{dataset}.csv"
        )

        performance.to_csv(
            performance_file,
            index=False
        )


    # --------------------------------------------------------
    # 8. Write readable analysis summary.
    # --------------------------------------------------------

    summary_file = (
        args.output_dir
        / f"analysis_{dataset}_summary.txt"
    )

    method_calls = (
        results
        .groupby("Method")["Significant"]
        .sum()
        .sort_index()
    )

    lines = [
        "DIFFERENTIAL-ABUNDANCE ANALYSIS SUMMARY",
        "=======================================",
        "",
        f"Dataset: {dataset}",
        f"Methods: {results['Method'].nunique()}",
        f"Unique taxa: {results['Taxon'].nunique()}",
        "",
        "SIGNIFICANT TAXA BY METHOD",
        "--------------------------"
    ]

    for method, count in method_calls.items():

        lines.append(
            f"{method}: {int(count)}"
        )

    lines.extend([
        "",
        "TAXON CONSENSUS",
        "---------------",
        (
            "Taxa called significant by at least one method: "
            f"{int((consensus['SignificantMethodCount'] > 0).sum())}"
        ),
        (
            "Taxa called significant by all methods: "
            f"{int((consensus['SignificantMethodCount'] == results['Method'].nunique()).sum())}"
        )
    ])


    # Ground-truth performance is only meaningful for simulated
    # datasets where the true DA status is known.
    if truth_available:

        lines.extend([
            "",
            "GROUND-TRUTH BENCHMARK",
            "----------------------",
            f"True DA taxa: {int(results.drop_duplicates('Taxon')['True_DA'].sum())}",
            ""
        ])

        for _, row in performance.iterrows():

            lines.extend([
                row["Method"],
                f"  True positives: {row['TruePositive']}",
                f"  False positives: {row['FalsePositive']}",
                f"  True negatives: {row['TrueNegative']}",
                f"  False negatives: {row['FalseNegative']}",
                f"  Specificity: {row['Specificity']:.4f}"
                if not pd.isna(row["Specificity"])
                else "  Specificity: NA",
                f"  False positive rate: {row['FalsePositiveRate']:.4f}"
                if not pd.isna(row["FalsePositiveRate"])
                else "  False positive rate: NA",
                ""
            ])


    summary_file.write_text(
        "\n".join(lines) + "\n"
    )


    # --------------------------------------------------------
    # 9. Terminal summary.
    # --------------------------------------------------------

    print()
    print("ANALYSIS COMPLETED")
    print("==================")

    print()
    print("Significant taxa by method:")

    for method, count in method_calls.items():
        print(
            f"  {method}: {int(count)}"
        )

    print()

    if truth_available:
        print(
            "Ground-truth performance table:"
        )
        print(
            f"  {performance_file}"
        )

    print(
        "Taxon consensus table:"
    )
    print(
        f"  {consensus_file}"
    )

    print(
        "Pairwise overlap table:"
    )
    print(
        f"  {overlap_file}"
    )

    print(
        "Analysis summary:"
    )
    print(
        f"  {summary_file}"
    )

    print()


# ============================================================
# RUN SCRIPT
# ============================================================

if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print()
        print("ANALYSIS FAILED")
        print("===============")
        print(error)
        print()

        sys.exit(1)
