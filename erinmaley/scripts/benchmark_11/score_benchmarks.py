from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("erinmaley/data/benchmark_11")
METHODS = ["ALDEx2", "ANCOM_BC", "CLR_Wilcoxon", "DESeq2"]
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
        results["Significant"].astype(str).str.lower() == "true"
    )

    for method in METHODS:
        subset = results[results["Method"] == method]

        calls = subset.set_index("Taxon")["Significant"]
        calls = calls.reindex(truth.index, fill_value=False)

        actual = truth["Planted_Target"].astype(bool)

        tp = int((calls & actual).sum())
        fp = int((calls & ~actual).sum())
        tn = int((~calls & ~actual).sum())
        fn = int((~calls & actual).sum())

        rows.append({
            "Simulation": i,
            "Method": method,
            "TP": tp,
            "FP": fp,
            "TN": tn,
            "FN": fn,
            "Sensitivity": tp / (tp + fn),
            "FDP": fp / (tp + fp) if tp + fp else np.nan,
            "FPR": fp / (fp + tn),
            "TaxaTested": len(subset)
        })

all_runs = pd.DataFrame(rows)

summary = all_runs.groupby("Method").agg(
    MeanSensitivity=("Sensitivity", "mean"),
    SDSensitivity=("Sensitivity", "std"),
    MeanFDP=("FDP", "mean"),
    SDFDP=("FDP", "std"),
    MeanFPR=("FPR", "mean"),
    MeanTaxaTested=("TaxaTested", "mean")
).reset_index()

all_runs.to_csv(BASE / "scored_all_runs.csv", index=False)
summary.to_csv(BASE / "scored_summary.csv", index=False)

print(summary.to_string(index=False))
print("SCORING COMPLETE")
