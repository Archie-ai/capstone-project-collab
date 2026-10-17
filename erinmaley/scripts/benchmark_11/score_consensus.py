from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("erinmaley/data/benchmark_11")
rows = []

for i in range(1, 6):
    folder = BASE / f"simulation_{i:03d}"

    truth = pd.read_csv(
        folder / "ground_truth.csv"
    ).set_index("Taxon")

    results = pd.read_csv(
        folder / "standard_results/standardized_signal_all_methods.csv"
    )

    results["Significant"] = (
        results["Significant"].astype(str).str.lower().eq("true")
    )

    significant = results[results["Significant"]]

    directional_counts = (
        significant.groupby(["Taxon", "Direction"])["Method"]
        .nunique()
        .unstack(fill_value=0)
    )

    for direction in ["Higher_in_Case", "Higher_in_Control"]:
        if direction not in directional_counts:
            directional_counts[direction] = 0

    actual = truth["Planted_Target"].astype(bool)

    for threshold in [2, 3, 4]:
        predicted = (
            directional_counts[
                ["Higher_in_Case", "Higher_in_Control"]
            ].max(axis=1)
            .reindex(truth.index, fill_value=0)
            .ge(threshold)
        )

        tp = int((predicted & actual).sum())
        fp = int((predicted & ~actual).sum())
        fn = int((~predicted & actual).sum())
        tn = int((~predicted & ~actual).sum())

        rows.append({
            "Simulation": i,
            "Threshold": threshold,
            "TP": tp,
            "FP": fp,
            "TN": tn,
            "FN": fn,
            "Sensitivity": tp / (tp + fn),
            "FDP": fp / (tp + fp) if tp + fp else np.nan,
            "FPR": fp / (fp + tn)
        })

all_runs = pd.DataFrame(rows)

summary = all_runs.groupby("Threshold").agg(
    MeanSensitivity=("Sensitivity", "mean"),
    SDSensitivity=("Sensitivity", "std"),
    MeanFDP=("FDP", "mean"),
    SDFDP=("FDP", "std"),
    MeanFPR=("FPR", "mean")
).reset_index()

all_runs.to_csv(BASE / "consensus_scored_all_runs.csv", index=False)
summary.to_csv(BASE / "consensus_scored_summary.csv", index=False)

print(summary.to_string(index=False))
print("CONSENSUS SCORING COMPLETE")
