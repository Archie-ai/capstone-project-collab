#!/usr/bin/env python3

# ============================================================
# CONSENSUS / REPORT
# ============================================================
#
# This script creates the final consensus report from the
# validated cross-method analysis.
#
# No differential-abundance testing is performed here.
#
# This script was developed with assistance from OpenAI's ChatGPT
#

import argparse
from pathlib import Path
import sys

import pandas as pd


# ============================================================
# COMMAND-LINE INPUTS
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="Create the final microbiome consensus report."
    )

    parser.add_argument(
        "--analysis-dir",
        required=True,
        type=Path,
        help="Stage 5 analysis output folder."
    )

    parser.add_argument(
        "--results-qc-dir",
        required=True,
        type=Path,
        help="Stage 6 results-QC folder."
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Folder for final consensus outputs."
    )

    parser.add_argument(
        "--dataset-label",
        required=True,
        help="Dataset label, such as null or real_test."
    )

    return parser.parse_args()


# ============================================================
# CONSENSUS CLASSIFICATION
# ============================================================

def classify_consensus(row):
    """
    Classify each taxon using the number of methods supporting it.

    Consensus means at least two methods called the taxon
    significant and those significant methods agree on direction.
    """

    count = int(
        row["SignificantMethodCount"]
    )

    direction_status = row[
        "DirectionAgreement"
    ]

    if count == 0:
        return "No_signal"

    if count == 1:
        return "Single_method_signal"

    if direction_status == "Agree":

        if count == 4:
            return "Full_consensus"

        if count == 3:
            return "Strong_consensus"

        return "Consensus"

    return "Discordant_signal"


# ============================================================
# CONSENSUS DIRECTION
# ============================================================

def get_consensus_direction(row):
    """Return a direction only when multiple methods agree."""

    if (
        row["SignificantMethodCount"] >= 2
        and row["DirectionAgreement"] == "Agree"
    ):

        return row["SignificantDirections"]

    return "Not_established"


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    dataset = args.dataset_label

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    print()
    print(
        f"Creating consensus report for dataset: {dataset}"
    )
    print()


    # --------------------------------------------------------
    # 1. Confirm Stage 6 Results QC passed.
    # --------------------------------------------------------

    qc_checks_file = (
        args.results_qc_dir
        / f"results_qc_{dataset}_checks.csv"
    )

    if not qc_checks_file.exists():
        raise FileNotFoundError(
            f"Results-QC checks not found: {qc_checks_file}"
        )

    qc_checks = pd.read_csv(
        qc_checks_file
    )

    failed_checks = qc_checks[
        qc_checks["Status"] != "PASS"
    ]

    if not failed_checks.empty:

        raise ValueError(
            "Stage 6 Results QC did not pass. "
            "Consensus report will not be created."
        )

    print(
        f"Results QC: PASS "
        f"({len(qc_checks)} checks)"
    )


    # --------------------------------------------------------
    # 2. Load Stage 5 taxon consensus table.
    # --------------------------------------------------------

    consensus_file = (
        args.analysis_dir
        / f"taxon_consensus_{dataset}.csv"
    )

    if not consensus_file.exists():
        raise FileNotFoundError(
            f"Taxon consensus table not found: {consensus_file}"
        )

    consensus = pd.read_csv(
        consensus_file
    )

    required_columns = {
        "Taxon",
        "MethodsSignificant",
        "SignificantMethodCount",
        "TotalMethods",
        "ConsensusFraction",
        "SignificantDirections",
        "DirectionAgreement"
    }

    missing = (
        required_columns
        - set(consensus.columns)
    )

    if missing:
        raise ValueError(
            "Consensus table is missing required column(s): "
            + ", ".join(sorted(missing))
        )


    # --------------------------------------------------------
    # 3. Assign the final consensus classification.
    # --------------------------------------------------------

    consensus["ConsensusClass"] = (
        consensus.apply(
            classify_consensus,
            axis=1
        )
    )

    consensus["ConsensusSignal"] = (
        consensus["ConsensusClass"].isin(
            [
                "Consensus",
                "Strong_consensus",
                "Full_consensus"
            ]
        )
    )

    consensus["ConsensusDirection"] = (
        consensus.apply(
            get_consensus_direction,
            axis=1
        )
    )


    # --------------------------------------------------------
    # 4. Create a clean final report table.
    # --------------------------------------------------------

    report_columns = [
        "Taxon",
        "MethodsSignificant",
        "SignificantMethodCount",
        "TotalMethods",
        "ConsensusFraction",
        "SignificantDirections",
        "DirectionAgreement",
        "ConsensusClass",
        "ConsensusSignal",
        "ConsensusDirection"
    ]

    # Simulation truth is retained when available.
    if "True_DA" in consensus.columns:
        report_columns.append(
            "True_DA"
        )

    if "True_Direction" in consensus.columns:
        report_columns.append(
            "True_Direction"
        )

    final_report = consensus[
        report_columns
    ].copy()


    # Put the strongest-supported taxa first.
    final_report = final_report.sort_values(
        by=[
            "SignificantMethodCount",
            "Taxon"
        ],
        ascending=[
            False,
            True
        ]
    )


    # --------------------------------------------------------
    # 5. Save full consensus report.
    # --------------------------------------------------------

    report_file = (
        args.output_dir
        / f"consensus_report_{dataset}.csv"
    )

    final_report.to_csv(
        report_file,
        index=False
    )


    # --------------------------------------------------------
    # 6. Save only consensus-positive taxa.
    # --------------------------------------------------------

    consensus_signals = final_report[
        final_report["ConsensusSignal"]
    ].copy()

    signal_file = (
        args.output_dir
        / f"consensus_signals_{dataset}.csv"
    )

    consensus_signals.to_csv(
        signal_file,
        index=False
    )


    # --------------------------------------------------------
    # 7. Create readable report summary.
    # --------------------------------------------------------

    class_counts = (
        final_report["ConsensusClass"]
        .value_counts()
    )

    lines = [
        "MICROBIOME CONSENSUS REPORT",
        "===========================",
        "",
        f"Dataset: {dataset}",
        "",
        "Results QC:",
        f"PASS ({len(qc_checks)} checks)",
        "",
        f"Taxa evaluated: {len(final_report)}",
        (
            "Taxa with at least one significant method: "
            f"{int((final_report['SignificantMethodCount'] > 0).sum())}"
        ),
        (
            "Consensus-positive taxa: "
            f"{int(final_report['ConsensusSignal'].sum())}"
        ),
        "",
        "CONSENSUS RULE",
        "--------------",
        (
            "A consensus signal requires at least two methods "
            "to call the taxon significant and agree on direction."
        ),
        "",
        "CONSENSUS CLASSES",
        "-----------------"
    ]

    class_order = [
        "Full_consensus",
        "Strong_consensus",
        "Consensus",
        "Discordant_signal",
        "Single_method_signal",
        "No_signal"
    ]

    for consensus_class in class_order:

        count = int(
            class_counts.get(
                consensus_class,
                0
            )
        )

        lines.append(
            f"{consensus_class}: {count}"
        )


    # --------------------------------------------------------
    # 8. Add benchmark information when simulation truth exists.
    # --------------------------------------------------------

    if "True_DA" in final_report.columns:

        true_da = (
            final_report["True_DA"]
            .astype(str)
            .str.lower()
            .map({
                "true": True,
                "false": False
            })
        )

        if true_da.notna().all():

            predicted = (
                final_report["ConsensusSignal"]
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

            lines.extend([
                "",
                "GROUND-TRUTH CONSENSUS BENCHMARK",
                "--------------------------------",
                f"True positives: {tp}",
                f"False positives: {fp}",
                f"True negatives: {tn}",
                f"False negatives: {fn}"
            ])


    # --------------------------------------------------------
    # 9. List consensus-positive taxa.
    # --------------------------------------------------------

    lines.extend([
        "",
        "CONSENSUS-POSITIVE TAXA",
        "-----------------------"
    ])

    if consensus_signals.empty:

        lines.append(
            "None"
        )

    else:

        for _, row in consensus_signals.iterrows():

            lines.append(
                (
                    f"{row['Taxon']}: "
                    f"{row['SignificantMethodCount']}/"
                    f"{row['TotalMethods']} methods; "
                    f"{row['ConsensusDirection']}; "
                    f"{row['MethodsSignificant']}"
                )
            )


    summary_file = (
        args.output_dir
        / f"consensus_report_{dataset}_summary.txt"
    )

    summary_file.write_text(
        "\n".join(lines) + "\n"
    )


    # --------------------------------------------------------
    # 10. Terminal output.
    # --------------------------------------------------------

    print()
    print("CONSENSUS REPORT COMPLETED")
    print("==========================")

    print(
        f"Taxa evaluated: {len(final_report)}"
    )

    print(
        "Consensus-positive taxa: "
        f"{int(final_report['ConsensusSignal'].sum())}"
    )

    print()
    print("Full report:")
    print(f"  {report_file}")

    print()
    print("Consensus signals:")
    print(f"  {signal_file}")

    print()
    print("Summary:")
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
        print("CONSENSUS REPORT FAILED")
        print("=======================")
        print(error)
        print()

        sys.exit(1)
