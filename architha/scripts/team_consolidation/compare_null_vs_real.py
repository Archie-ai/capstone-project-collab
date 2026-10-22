#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUT_DIR = Path(
    "architha/results/team_consolidation/null_vs_real"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def pct(numerator, denominator):
    if numerator is None or denominator is None:
        return None

    if denominator == 0:
        return None

    return 100.0 * numerator / denominator


def add_method(
    rows,
    dataset,
    owner,
    method,
    variant,
    universe,
    valid_q,
    significant,
    robust=None,
    notes=""
):
    rows.append(
        {
            "Dataset": dataset,
            "Owner": owner,
            "Method": method,
            "Method_Variant": variant,
            "Analysis_Universe": universe,
            "Valid_Adjusted_P_Values": valid_q,
            "Significant_Taxa": significant,
            "Discovery_Rate_Universe_Pct": pct(
                significant,
                universe
            ),
            "Discovery_Rate_ValidQ_Pct": pct(
                significant,
                valid_q
            ),
            "Robust_Significant_Taxa": robust,
            "Notes": notes,
        }
    )


# ============================================================
# STORAGE
# ============================================================

method_rows = []

dataset_rows = []


# ============================================================
# 1. CANONICAL NULL
# Verified independent null validation
# 50 samples, 100 taxa, 0 true DA taxa
# ============================================================

add_method(
    method_rows,
    dataset="canonical_null",
    owner="Architha",
    method="CLR_Wilcoxon",
    variant="CLR + Wilcoxon/Mann-Whitney",
    universe=100,
    valid_q=100,
    significant=1,
    notes="Known complete-null dataset."
)

add_method(
    method_rows,
    dataset="canonical_null",
    owner="Architha",
    method="ALDEx2",
    variant="ALDEx2 Wilcoxon primary",
    universe=100,
    valid_q=100,
    significant=0,
    notes="Known complete-null dataset."
)

add_method(
    method_rows,
    dataset="canonical_null",
    owner="Architha",
    method="ANCOM",
    variant="ANCOM-BC2 robust primary",
    universe=100,
    valid_q=100,
    significant=1,
    robust=1,
    notes="Known complete-null dataset."
)

add_method(
    method_rows,
    dataset="canonical_null",
    owner="Architha",
    method="DESeq2",
    variant="DESeq2 primary",
    universe=100,
    valid_q=100,
    significant=0,
    notes="Known complete-null dataset."
)

dataset_rows.append(
    {
        "Dataset": "canonical_null",
        "Owner": "Architha",
        "Dataset_Type": "Known complete null",
        "Samples": 50,
        "Group1_N": 25,
        "Group2_N": 25,
        "Raw_Taxa": 100,
        "QC_Retained_Taxa": 100,
        "Analysis_Universe": 100,
        "Zero_Fraction": 0.0468,
        "Preferred_Consensus_ge2": 1,
        "Preferred_Consensus_ge3": 0,
        "Preferred_Consensus_4of4": 0,
        "Consensus_Definition": "Vote-based null consensus",
        "ANCOM_Standard_Significant": 1,
        "ANCOM_Robust_Significant": 1,
        "Notes": (
            "Ground truth contains zero true DA taxa. "
            "Taxon_033 was a 2-of-4 false consensus."
        ),
    }
)


# ============================================================
# 2. ERIN - CDI_SCHUBERT
# Verified current repository results
# ============================================================

add_method(
    method_rows,
    dataset="cdi_schubert",
    owner="Erin",
    method="CLR_Wilcoxon",
    variant="CLR + Wilcoxon; pseudocount 1.0",
    universe=1420,
    valid_q=1420,
    significant=1374
)

add_method(
    method_rows,
    dataset="cdi_schubert",
    owner="Erin",
    method="ALDEx2",
    variant="ALDEx2 Wilcoxon primary; 128 MC",
    universe=1420,
    valid_q=1420,
    significant=92
)

add_method(
    method_rows,
    dataset="cdi_schubert",
    owner="Erin",
    method="ANCOM",
    variant="ANCOM-BC standard",
    universe=1420,
    valid_q=1420,
    significant=521,
    notes="Robust ANCOM sensitivity not reported."
)

add_method(
    method_rows,
    dataset="cdi_schubert",
    owner="Erin",
    method="DESeq2",
    variant="DESeq2 primary",
    universe=1420,
    valid_q=457,
    significant=283,
    notes="963 taxa had NA adjusted p-values."
)

dataset_rows.append(
    {
        "Dataset": "cdi_schubert",
        "Owner": "Erin",
        "Dataset_Type": "Real",
        "Samples": 247,
        "Group1_N": 154,
        "Group2_N": 93,
        "Raw_Taxa": 19314,
        "QC_Retained_Taxa": 1420,
        "Analysis_Universe": 1420,
        "Zero_Fraction": 0.8690,
        "Preferred_Consensus_ge2": 286,
        "Preferred_Consensus_ge3": 106,
        "Preferred_Consensus_4of4": 59,
        "Consensus_Definition": "Direction-aware consensus",
        "ANCOM_Standard_Significant": 521,
        "ANCOM_Robust_Significant": None,
        "Notes": (
            "Preferred consensus requires significance "
            "and agreement in effect direction."
        ),
    }
)


# ============================================================
# 3. PEGGY - CRC_BAXTER
# Method results verified from repository.
# Vote consensus calculated from Peggy's final consensus file.
# ============================================================

peggy_metadata = pd.read_csv(
    "peggyhall/data/real/crc_baxter/"
    "crc_baxter_filtered_metadata.csv"
)

peggy_consensus = pd.read_csv(
    "peggyhall/data/real/crc_baxter/DA_results/"
    "crc_baxter_consensus_all_taxa.csv"
)

peggy_votes = pd.to_numeric(
    peggy_consensus["Votes"],
    errors="coerce"
)

peggy_ge2 = int(
    (peggy_votes >= 2).sum()
)

peggy_ge3 = int(
    (peggy_votes >= 3).sum()
)

peggy_4 = int(
    (peggy_votes == 4).sum()
)

add_method(
    method_rows,
    dataset="crc_baxter",
    owner="Peggy",
    method="CLR_Wilcoxon",
    variant="CLR + Wilcoxon",
    universe=1133,
    valid_q=1133,
    significant=24
)

add_method(
    method_rows,
    dataset="crc_baxter",
    owner="Peggy",
    method="ALDEx2",
    variant="ALDEx2 Welch primary",
    universe=1133,
    valid_q=1133,
    significant=0,
    notes=(
        "Welch-primary ALDEx2 differs from "
        "the Wilcoxon-primary workflows."
    )
)

add_method(
    method_rows,
    dataset="crc_baxter",
    owner="Peggy",
    method="ANCOM",
    variant="ANCOM-BC2 standard",
    universe=1133,
    valid_q=1133,
    significant=486,
    notes="Robust ANCOM sensitivity not reported."
)

add_method(
    method_rows,
    dataset="crc_baxter",
    owner="Peggy",
    method="DESeq2",
    variant="DESeq2 primary",
    universe=1133,
    valid_q=584,
    significant=20,
    notes="549 taxa lacked usable adjusted p-values."
)

dataset_rows.append(
    {
        "Dataset": "crc_baxter",
        "Owner": "Peggy",
        "Dataset_Type": "Real",
        "Samples": len(peggy_metadata),
        "Group1_N": None,
        "Group2_N": None,
        "Raw_Taxa": None,
        "QC_Retained_Taxa": 1133,
        "Analysis_Universe": 1133,
        "Zero_Fraction": None,
        "Preferred_Consensus_ge2": peggy_ge2,
        "Preferred_Consensus_ge3": peggy_ge3,
        "Preferred_Consensus_4of4": peggy_4,
        "Consensus_Definition": "Vote-based consensus",
        "ANCOM_Standard_Significant": 486,
        "ANCOM_Robust_Significant": None,
        "Notes": "ALDEx2 is Welch-primary.",
    }
)


# ============================================================
# 4. RAQUEL - CDI_VINCENT
# Verified current repository results
# ============================================================

add_method(
    method_rows,
    dataset="cdi_vincent",
    owner="Raquel",
    method="CLR_Wilcoxon",
    variant="CLR + Wilcoxon; pseudocount 0.5",
    universe=195,
    valid_q=195,
    significant=14
)

add_method(
    method_rows,
    dataset="cdi_vincent",
    owner="Raquel",
    method="ALDEx2",
    variant="ALDEx2 Wilcoxon primary; 128 MC",
    universe=195,
    valid_q=195,
    significant=3
)

add_method(
    method_rows,
    dataset="cdi_vincent",
    owner="Raquel",
    method="ANCOM",
    variant="ANCOM-BC2 standard + robust",
    universe=195,
    valid_q=195,
    significant=61,
    robust=1,
    notes="Standard and robust ANCOM calls retained."
)

add_method(
    method_rows,
    dataset="cdi_vincent",
    owner="Raquel",
    method="DESeq2",
    variant="DESeq2 primary",
    universe=195,
    valid_q=142,
    significant=34,
    notes="53 taxa had missing adjusted p-values."
)

dataset_rows.append(
    {
        "Dataset": "cdi_vincent",
        "Owner": "Raquel",
        "Dataset_Type": "Real",
        "Samples": 50,
        "Group1_N": 25,
        "Group2_N": 25,
        "Raw_Taxa": 4299,
        "QC_Retained_Taxa": 618,
        "Analysis_Universe": 195,
        "Zero_Fraction": 0.806990291262136,
        "Preferred_Consensus_ge2": 6,
        "Preferred_Consensus_ge3": 3,
        "Preferred_Consensus_4of4": 1,
        "Consensus_Definition": "Robust vote-based consensus",
        "ANCOM_Standard_Significant": 61,
        "ANCOM_Robust_Significant": 1,
        "Notes": (
            "Preprocessing retained 618 OTUs; "
            "final DA summaries use 195 OTUs."
        ),
    }
)


# ============================================================
# 5. ARCHITHA - PAR_SCHEPERJANS
# Verified current repository results
# ============================================================

add_method(
    method_rows,
    dataset="par_scheperjans",
    owner="Architha",
    method="CLR_Wilcoxon",
    variant=(
        "CLR + Mann-Whitney/Wilcoxon; pseudocount 0.5"
    ),
    universe=691,
    valid_q=691,
    significant=1
)

add_method(
    method_rows,
    dataset="par_scheperjans",
    owner="Architha",
    method="ALDEx2",
    variant="ALDEx2 Wilcoxon primary; 128 MC",
    universe=691,
    valid_q=691,
    significant=0
)

add_method(
    method_rows,
    dataset="par_scheperjans",
    owner="Architha",
    method="ANCOM",
    variant="ANCOM-BC2 standard + robust",
    universe=691,
    valid_q=691,
    significant=365,
    robust=1,
    notes=(
        "Standard calls used for harmonized comparison; "
        "robust sensitivity retained separately."
    )
)

add_method(
    method_rows,
    dataset="par_scheperjans",
    owner="Architha",
    method="DESeq2",
    variant="DESeq2 poscounts; local dispersion fit",
    universe=691,
    valid_q=691,
    significant=3,
    notes=(
        "Primary and diagnostic DESeq2 analyses "
        "all returned 3 significant taxa."
    )
)

dataset_rows.append(
    {
        "Dataset": "par_scheperjans",
        "Owner": "Architha",
        "Dataset_Type": "Real",
        "Samples": 148,
        "Group1_N": 74,
        "Group2_N": 74,
        "Raw_Taxa": 10232,
        "QC_Retained_Taxa": 691,
        "Analysis_Universe": 691,
        "Zero_Fraction": 0.784458,
        "Preferred_Consensus_ge2": 1,
        "Preferred_Consensus_ge3": 1,
        "Preferred_Consensus_4of4": 0,
        "Consensus_Definition": (
            "Robust ANCOM + direction-aware consensus"
        ),
        "ANCOM_Standard_Significant": 365,
        "ANCOM_Robust_Significant": 1,
        "Notes": (
            "OTU_00659 / Alistipes is the "
            "single 3-of-4 consensus taxon."
        ),
    }
)


# ============================================================
# 6. CURRENT 200-RUN NULL CALIBRATION
# Verified from current shared repository
# ============================================================

null_method_rows = [
    {
        "Method": "CLR_Wilcoxon",
        "Mean_Type_I_Error": 0.00035,
        "Probability_Any_False_Discovery": 0.035,
        "Probability_Any_False_Discovery_Pct": 3.5,
        "Mean_False_Positives": 0.035,
    },
    {
        "Method": "ALDEx2",
        "Mean_Type_I_Error": 0.00030,
        "Probability_Any_False_Discovery": 0.030,
        "Probability_Any_False_Discovery_Pct": 3.0,
        "Mean_False_Positives": 0.030,
    },
    {
        "Method": "ANCOM",
        "Mean_Type_I_Error": 0.00135,
        "Probability_Any_False_Discovery": 0.095,
        "Probability_Any_False_Discovery_Pct": 9.5,
        "Mean_False_Positives": 0.135,
    },
    {
        "Method": "DESeq2",
        "Mean_Type_I_Error": 0.00180,
        "Probability_Any_False_Discovery": 0.140,
        "Probability_Any_False_Discovery_Pct": 14.0,
        "Mean_False_Positives": 0.180,
    },
]

null_consensus_rows = [
    {
        "Threshold": "2-of-4",
        "Mean_Type_I_Error": 0.0008,
        "Probability_Any_False_Discovery": 0.07,
        "Probability_Any_False_Discovery_Pct": 7.0,
        "Mean_False_Positives": 0.08,
    },
    {
        "Threshold": "3-of-4",
        "Mean_Type_I_Error": 0.0003,
        "Probability_Any_False_Discovery": 0.03,
        "Probability_Any_False_Discovery_Pct": 3.0,
        "Mean_False_Positives": 0.03,
    },
    {
        "Threshold": "4-of-4",
        "Mean_Type_I_Error": 0.0002,
        "Probability_Any_False_Discovery": 0.02,
        "Probability_Any_False_Discovery_Pct": 2.0,
        "Mean_False_Positives": 0.02,
    },
]


# ============================================================
# 7. BUILD TABLES
# ============================================================

method_df = pd.DataFrame(
    method_rows
)

dataset_df = pd.DataFrame(
    dataset_rows
)

null_method_df = pd.DataFrame(
    null_method_rows
)

null_consensus_df = pd.DataFrame(
    null_consensus_rows
)


dataset_order = [
    "canonical_null",
    "cdi_schubert",
    "crc_baxter",
    "cdi_vincent",
    "par_scheperjans",
]

method_order = [
    "CLR_Wilcoxon",
    "ALDEx2",
    "ANCOM",
    "DESeq2",
]


method_df["Dataset"] = pd.Categorical(
    method_df["Dataset"],
    categories=dataset_order,
    ordered=True
)

method_df["Method"] = pd.Categorical(
    method_df["Method"],
    categories=method_order,
    ordered=True
)

method_df = (
    method_df
    .sort_values(
        ["Dataset", "Method"]
    )
    .reset_index(drop=True)
)


dataset_df["Dataset"] = pd.Categorical(
    dataset_df["Dataset"],
    categories=dataset_order,
    ordered=True
)

dataset_df = (
    dataset_df
    .sort_values("Dataset")
    .reset_index(drop=True)
)


# ============================================================
# 8. VALIDATION
# ============================================================

if len(method_df) != 20:
    raise ValueError(
        f"Expected 20 method rows, found {len(method_df)}."
    )

if len(dataset_df) != 5:
    raise ValueError(
        f"Expected 5 dataset rows, found {len(dataset_df)}."
    )


# ============================================================
# 9. WIDE DISCOVERY-RATE TABLE
# ============================================================

wide_df = (
    method_df
    .pivot(
        index="Dataset",
        columns="Method",
        values="Discovery_Rate_Universe_Pct"
    )
    .reset_index()
)


# ============================================================
# 10. SAVE OUTPUTS
# ============================================================

method_df.to_csv(
    OUT_DIR / "null_vs_real_method_comparison.csv",
    index=False
)

dataset_df.to_csv(
    OUT_DIR / "null_vs_real_dataset_summary.csv",
    index=False
)

wide_df.to_csv(
    OUT_DIR / "null_vs_real_discovery_rates_wide.csv",
    index=False
)

null_method_df.to_csv(
    OUT_DIR / "null_200_method_calibration_current.csv",
    index=False
)

null_consensus_df.to_csv(
    OUT_DIR / "null_200_consensus_calibration_current.csv",
    index=False
)


# ============================================================
# 11. INTERPRETATION FILE
# ============================================================

report = [
    "NULL-VERSUS-REAL MICROBIOME DA CONSOLIDATION",
    "=" * 68,
    "",
    "The canonical null simulation contains zero true DA taxa.",
    (
        "Therefore, any discovery in the null dataset is a "
        "false discovery for that realization."
    ),
    "",
    (
        "The real microbiome datasets do not have known "
        "biological ground truth."
    ),
    (
        "Real-data significant taxa and consensus signals "
        "must therefore be interpreted as statistical "
        "discoveries and cross-method agreement, not proven "
        "true positives."
    ),
    "",
    "200-RUN NULL CALIBRATION",
    "------------------------",
    (
        "CLR/Wilcoxon: 3.5% of runs contained at least "
        "one false discovery."
    ),
    (
        "ALDEx2: 3.0% of runs contained at least "
        "one false discovery."
    ),
    (
        "ANCOM-BC: 9.5% of runs contained at least "
        "one false discovery."
    ),
    (
        "DESeq2: 14.0% of runs contained at least "
        "one false discovery."
    ),
    "",
    "CONSENSUS NULL CALIBRATION",
    "--------------------------",
    (
        "2-of-4 consensus: 7.0% probability of at least "
        "one false consensus discovery."
    ),
    (
        "3-of-4 consensus: 3.0% probability of at least "
        "one false consensus discovery."
    ),
    (
        "4-of-4 consensus: 2.0% probability of at least "
        "one false consensus discovery."
    ),
    "",
    "COMPARABILITY NOTES",
    "-------------------",
    (
        "Discovery rates are normalized by each dataset's "
        "analysis universe because raw significant-taxon "
        "counts are not directly comparable."
    ),
    (
        "DESeq2 valid-q discovery rates are retained "
        "because some datasets contain many NA adjusted "
        "p-values."
    ),
    (
        "Standard ANCOM calls are used for cross-dataset "
        "comparison. Robust ANCOM calls are retained "
        "separately where available."
    ),
    (
        "Peggy's ALDEx2 analysis is Welch-primary, "
        "whereas the other ALDEx2 workflows are "
        "Wilcoxon-primary."
    ),
    (
        "Consensus definitions differ somewhat among "
        "team workflows and are explicitly labeled."
    ),
    "",
    "MAIN INTERPRETATION",
    "-------------------",
    (
        "Discovery burden varies strongly across both "
        "methods and real datasets."
    ),
    (
        "The magnitude of this variation is substantially "
        "greater than the average false-positive burden "
        "observed under repeated complete-null simulation."
    ),
    (
        "However, because biological truth is unknown in "
        "the real datasets, larger discovery sets cannot "
        "automatically be interpreted as greater accuracy."
    ),
    (
        "The results instead demonstrate strong "
        "method-by-dataset dependence."
    ),
    (
        "Consensus provides a way to prioritize more "
        "reproducible signals, but the null results show "
        "that consensus itself is not equivalent to truth."
    ),
]


(
    OUT_DIR / "null_vs_real_interpretation.txt"
).write_text(
    "\n".join(report),
    encoding="utf-8"
)


# ============================================================
# 12. TERMINAL OUTPUT
# ============================================================

print("=" * 72)
print("NULL-VERSUS-REAL TEAM CONSOLIDATION")
print("=" * 72)

print("\nDATASET SUMMARY\n")

print(
    dataset_df[
        [
            "Dataset",
            "Owner",
            "Samples",
            "Analysis_Universe",
            "Preferred_Consensus_ge2",
            "Preferred_Consensus_ge3",
            "Preferred_Consensus_4of4",
        ]
    ].to_string(index=False)
)

print("\nMETHOD COMPARISON\n")

print(
    method_df[
        [
            "Dataset",
            "Method",
            "Analysis_Universe",
            "Valid_Adjusted_P_Values",
            "Significant_Taxa",
            "Discovery_Rate_Universe_Pct",
            "Discovery_Rate_ValidQ_Pct",
            "Robust_Significant_Taxa",
        ]
    ].to_string(index=False)
)

print("\n200-RUN NULL METHOD CALIBRATION\n")

print(
    null_method_df.to_string(
        index=False
    )
)

print("\n200-RUN NULL CONSENSUS CALIBRATION\n")

print(
    null_consensus_df.to_string(
        index=False
    )
)

print("\nCreated outputs:")

for filename in [
    "null_vs_real_method_comparison.csv",
    "null_vs_real_dataset_summary.csv",
    "null_vs_real_discovery_rates_wide.csv",
    "null_200_method_calibration_current.csv",
    "null_200_consensus_calibration_current.csv",
    "null_vs_real_interpretation.txt",
]:
    print(
        OUT_DIR / filename
    )

print("\n" + "=" * 72)
print("NULL-VERSUS-REAL TEAM CONSOLIDATION: PASS")
print("=" * 72)
