#!/usr/bin/env python3

# ============================================================
# STANDARD RESULTS TABLES
# ============================================================
#
# This script converts the four differential-abundance method
# outputs into one shared format for downstream analysis.
#
# Original method outputs are not changed.
#
# This script was developed with assistance from OpenAI's ChatGPT
#

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd


# ============================================================
# COMMAND-LINE INPUTS
# ============================================================

def parse_args():
    """Collect the method-output folder, output folder, and dataset name."""

    parser = argparse.ArgumentParser(
        description="Standardize results from all four DA methods."
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Folder containing the four method output folders."
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Folder where standardized results will be saved."
    )

    parser.add_argument(
        "--dataset-label",
        required=True,
        help="Dataset label, such as null or real_test."
    )

    return parser.parse_args()


# ============================================================
# STANDARDIZE ONE METHOD
# ============================================================

def standardize_method(
    data,
    dataset,
    method,
    effect_column,
    effect_metric,
    statistic_column=None
):
    """
    Convert one method's output into the shared result format.

    Effect sizes are preserved, but EffectMetric records what each
    number actually represents because effect-size scales differ
    between methods.
    """

    required = {
        "Taxon",
        effect_column,
        "PValue",
        "AdjustedPValue",
        "Significant",
        "Direction"
    }

    missing = required - set(data.columns)

    if missing:
        raise ValueError(
            f"{method} is missing required column(s): "
            + ", ".join(sorted(missing))
        )

    if statistic_column is not None:

        if statistic_column not in data.columns:
            raise ValueError(
                f"{method} is missing statistic column: "
                f"{statistic_column}"
            )

        statistic = data[statistic_column]

    else:
        # ALDEx2 does not have the same single test-statistic
        # field used by the other methods.
        statistic = np.nan


    standardized = pd.DataFrame({
        "Dataset": dataset,
        "Method": method,
        "Taxon": data["Taxon"],
        "EffectSize": data[effect_column],
        "EffectMetric": effect_metric,
        "Statistic": statistic,
        "PValue": data["PValue"],
        "AdjustedPValue": data["AdjustedPValue"],
        "Significant": data["Significant"],
        "Direction": data["Direction"]
    })

    return standardized


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    dataset = args.dataset_label

    print()
    print(f"Standardizing results for dataset: {dataset}")
    print()


    # --------------------------------------------------------
    # 1. Locate the four Stage 3 result files.
    # --------------------------------------------------------

    clr_file = (
        args.input_dir
        / "clr_wilcoxon"
        / f"clr_wilcoxon_{dataset}_results.csv"
    )

    aldex2_file = (
        args.input_dir
        / "aldex2"
        / f"aldex2_{dataset}_results.csv"
    )

    ancombc_file = (
        args.input_dir
        / "ancombc"
        / f"ancombc_{dataset}_results.csv"
    )

    deseq2_file = (
        args.input_dir
        / "deseq2"
        / f"deseq2_{dataset}_results.csv"
    )

    input_files = [
        clr_file,
        aldex2_file,
        ancombc_file,
        deseq2_file
    ]

    for file in input_files:

        if not file.exists():
            raise FileNotFoundError(
                f"Method result file not found: {file}"
            )


    # --------------------------------------------------------
    # 2. Load the original method outputs.
    # --------------------------------------------------------

    clr = pd.read_csv(clr_file)
    aldex2 = pd.read_csv(aldex2_file)
    ancombc = pd.read_csv(ancombc_file)
    deseq2 = pd.read_csv(deseq2_file)


    # --------------------------------------------------------
    # 3. Standardize CLR/Wilcoxon.
    # --------------------------------------------------------

    clr_standard = standardize_method(
        data=clr,
        dataset=dataset,
        method="CLR_Wilcoxon",
        effect_column="CLRDifference",
        effect_metric="CLR_median_difference",
        statistic_column="Statistic"
    )


    # --------------------------------------------------------
    # 4. Standardize ALDEx2.
    # --------------------------------------------------------

    aldex2_standard = standardize_method(
        data=aldex2,
        dataset=dataset,
        method="ALDEx2",
        effect_column="DifferenceBetween",
        effect_metric="CLR_difference_between",
        statistic_column=None
    )


    # --------------------------------------------------------
    # 5. Standardize ANCOM-BC.
    # --------------------------------------------------------

    ancombc_standard = standardize_method(
        data=ancombc,
        dataset=dataset,
        method="ANCOM_BC",
        effect_column="LogFoldChange",
        effect_metric="log_fold_change",
        statistic_column="Statistic"
    )


    # --------------------------------------------------------
    # 6. Standardize DESeq2.
    # --------------------------------------------------------

    deseq2_standard = standardize_method(
        data=deseq2,
        dataset=dataset,
        method="DESeq2",
        effect_column="Log2FoldChange",
        effect_metric="log2_fold_change",
        statistic_column="Statistic"
    )


    # --------------------------------------------------------
    # 7. Save each standardized method table.
    # --------------------------------------------------------

    method_tables = {
        "clr_wilcoxon": clr_standard,
        "aldex2": aldex2_standard,
        "ancombc": ancombc_standard,
        "deseq2": deseq2_standard
    }

    for method_name, table in method_tables.items():

        output_file = (
            args.output_dir
            / f"standardized_{method_name}_{dataset}.csv"
        )

        table.to_csv(
            output_file,
            index=False
        )

        print(
            f"{method_name}: "
            f"{len(table)} taxa standardized"
        )


    # --------------------------------------------------------
    # 8. Combine all four methods into one long-format table.
    # --------------------------------------------------------

    combined = pd.concat(
        [
            clr_standard,
            aldex2_standard,
            ancombc_standard,
            deseq2_standard
        ],
        ignore_index=True
    )

    combined_file = (
        args.output_dir
        / f"standardized_{dataset}_all_methods.csv"
    )

    combined.to_csv(
        combined_file,
        index=False
    )


    # --------------------------------------------------------
    # 9. Confirm that all methods contributed results.
    # --------------------------------------------------------

    expected_methods = 4

    if combined["Method"].nunique() != expected_methods:

        raise ValueError(
            "Combined table does not contain all four methods."
        )


    # --------------------------------------------------------
    # 10. Terminal summary.
    # --------------------------------------------------------

    print()
    print("STANDARDIZATION COMPLETED")
    print("=========================")

    print(
        f"Methods: {combined['Method'].nunique()}"
    )

    print(
        f"Total rows: {len(combined)}"
    )

    print(
        f"Unique taxa: {combined['Taxon'].nunique()}"
    )

    print()
    print("Combined table:")
    print(f"  {combined_file}")

    print()


# ============================================================
# RUN SCRIPT
# ============================================================

if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print()
        print("STANDARDIZATION FAILED")
        print("======================")
        print(error)
        print()

        sys.exit(1)
