from pathlib import Path
import pandas as pd

BASE = Path("erinmaley/data/null_benchmark_12")
METHODS = ["ALDEx2", "ANCOM_BC", "CLR_Wilcoxon", "DESeq2"]

method_rows = []
consensus_rows = []

for i in range(1, 201):
    folder = BASE / f"simulation_{i:03d}"

    truth = pd.read_csv(folder / "ground_truth.csv")
    taxa = pd.Index(truth["Taxon"])

    assert len(taxa) == 100
    assert not truth["True_DA"].astype(bool).any()

    results = pd.read_csv(
        folder / "standard_results/standardized_signal_all_methods.csv"
    )

    results["Significant"] = (
        results["Significant"].astype(str).str.lower().eq("true")
    )

    assert set(results["Method"]) == set(METHODS)
    assert not results.duplicated(["Method", "Taxon"]).any()
    assert set(results["Taxon"]).issubset(set(taxa))

    # Individual method performance
    for method in METHODS:
        subset = results[results["Method"] == method]

        calls = subset.set_index("Taxon")["Significant"]
        calls = calls.reindex(taxa, fill_value=False)

        fp = int(calls.sum())

        method_rows.append({
            "Simulation": i,
            "Method": method,
            "FalsePositives": fp,
            "TypeIError": fp / len(taxa),
            "AnyFalseDiscovery": int(fp > 0)
        })

    # Direction-aware consensus voting
    significant = results[results["Significant"]]

    direction_counts = (
        significant.groupby(["Taxon", "Direction"])["Method"]
        .nunique()
        .unstack(fill_value=0)
    )

    for direction in ["Higher_in_Case", "Higher_in_Control"]:
        if direction not in direction_counts.columns:
            direction_counts[direction] = 0

    max_votes = (
        direction_counts[
            ["Higher_in_Case", "Higher_in_Control"]
        ]
        .max(axis=1)
        .reindex(taxa, fill_value=0)
    )

    for threshold in [2, 3, 4]:
        fp = int((max_votes >= threshold).sum())

        consensus_rows.append({
            "Simulation": i,
            "Threshold": threshold,
            "FalsePositives": fp,
            "TypeIError": fp / len(taxa),
            "AnyFalseDiscovery": int(fp > 0)
        })

methods = pd.DataFrame(method_rows)
consensus = pd.DataFrame(consensus_rows)

method_summary = methods.groupby("Method").agg(
    MeanTypeIError=("TypeIError", "mean"),
    ProbabilityAnyFalseDiscovery=("AnyFalseDiscovery", "mean"),
    MeanFalsePositives=("FalsePositives", "mean")
).reset_index()

consensus_summary = consensus.groupby("Threshold").agg(
    MeanTypeIError=("TypeIError", "mean"),
    ProbabilityAnyFalseDiscovery=("AnyFalseDiscovery", "mean"),
    MeanFalsePositives=("FalsePositives", "mean")
).reset_index()

methods.to_csv(BASE / "null_200_method_all_runs.csv", index=False)
method_summary.to_csv(BASE / "null_200_method_summary.csv", index=False)

consensus.to_csv(BASE / "null_200_consensus_all_runs.csv", index=False)
consensus_summary.to_csv(BASE / "null_200_consensus_summary.csv", index=False)

print("\nMETHOD PERFORMANCE — 200 SIMULATIONS")
print(method_summary.to_string(index=False))

print("\nCONSENSUS PERFORMANCE — 200 SIMULATIONS")
print(consensus_summary.to_string(index=False))

print("\n200-SIMULATION NULL BENCHMARK SCORING COMPLETE")
