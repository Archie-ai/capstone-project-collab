# BIOT 670I - Microbiome Abundance Consensus Machine
# Consensus analysis across 100 Version 2 simulations

from pathlib import Path
import numpy as np
import pandas as pd

N_SIMULATIONS = 100

BASE = Path("results/da")
OUTPUT = BASE / "consensus_v2"
OUTPUT.mkdir(parents=True, exist_ok=True)

METHODS = {
    "CLR_Wilcoxon": ("clr_wilcoxon_v2", "Significant"),
    "ALDEx2": ("aldex2_v2", "Significant"),
    "ANCOMBC2": ("ancombc_v2", "Robust_Significant"),
    "DESeq2": ("deseq2_v2", "Significant"),
}

RULES = [2, 3, 4]

performance_rows = []
agreement_rows = []

for sim in range(1, N_SIMULATIONS + 1):

    sim_name = f"simulation_{sim:03d}"

    combined = None

    for method, (folder, sig_col) in METHODS.items():

        path = BASE / folder / f"{sim_name}_results.csv"
        df = pd.read_csv(path)

        required = {"Taxon", "True_DA", sig_col}
        if not required.issubset(df.columns):
            raise ValueError(f"Missing columns in {path}")

        if df["Taxon"].duplicated().any():
            raise ValueError(f"Duplicate taxa in {path}")

        if df["True_DA"].isna().any():
            raise ValueError(f"Missing ground truth in {path}")

        # Missing significance means the taxon was not detected
        detected = (
            df[sig_col]
            .fillna(False)
            .astype(bool)
        )

        current = pd.DataFrame({
            "Taxon": df["Taxon"],
            "True_DA": df["True_DA"].astype(bool),
            method: detected,
        })

        if combined is None:
            combined = current
        else:
            combined = combined.merge(
                current,
                on=["Taxon", "True_DA"],
                how="outer",
                validate="one_to_one"
            )

    method_cols = list(METHODS)

    # Missing method calls count as not detected
    combined[method_cols] = (
        combined[method_cols].fillna(False).astype(bool)
    )

    combined["Votes"] = combined[method_cols].sum(axis=1)

    # Save method-level detections for future figures
    combined.to_csv(
        OUTPUT / f"{sim_name}_votes.csv",
        index=False
    )

    # Number of taxa receiving 0, 1, 2, 3, or 4 votes
    vote_counts = combined["Votes"].value_counts()

    for votes in range(5):
        agreement_rows.append({
            "Simulation": sim,
            "Votes": votes,
            "Taxa_Count": int(vote_counts.get(votes, 0)),
        })

    truth = combined["True_DA"].to_numpy()

    for minimum_votes in RULES:

        detected = (
            combined["Votes"].to_numpy() >= minimum_votes
        )

        TP = int(np.sum(detected & truth))
        FP = int(np.sum(detected & ~truth))
        FN = int(np.sum(~detected & truth))
        TN = int(np.sum(~detected & ~truth))

        discoveries = TP + FP

        performance_rows.append({
            "Simulation": sim,
            "Rule": f"{minimum_votes}-of-4",
            "TP": TP,
            "FP": FP,
            "FN": FN,
            "TN": TN,
            "Discoveries": discoveries,
            "FDP": FP / discoveries if discoveries else 0.0,
            "Power": TP / (TP + FN) if TP + FN else np.nan,
            "Sensitivity": TP / (TP + FN) if TP + FN else np.nan,
            "Specificity": TN / (TN + FP) if TN + FP else np.nan,
        })

    print(f"Simulation {sim}/{N_SIMULATIONS} analyzed.")

# Save per-simulation consensus performance
performance = pd.DataFrame(performance_rows)

performance.to_csv(
    OUTPUT / "consensus_performance_by_simulation.csv",
    index=False
)

# Summarize each voting rule
metrics = [
    "TP", "FP", "FN", "TN",
    "Discoveries", "FDP",
    "Power", "Sensitivity", "Specificity"
]

summary = (
    performance
    .groupby("Rule")[metrics]
    .agg(["mean", "std", "min", "max"])
)

summary.columns = [
    f"{metric}_{stat}"
    for metric, stat in summary.columns
]

summary = summary.reset_index()

summary.to_csv(
    OUTPUT / "consensus_performance_summary.csv",
    index=False
)

# Save voting agreement distribution
agreement = pd.DataFrame(agreement_rows)

agreement.to_csv(
    OUTPUT / "vote_distribution_by_simulation.csv",
    index=False
)

print("\nAll 100 consensus analyses completed!")
print(
    summary[
        ["Rule", "FDP_mean", "Power_mean", "Specificity_mean"]
    ].to_string(index=False)
)

print("\nResults saved to:", OUTPUT)
