import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

# Settings
n_simulations = 100
pseudocount = 0.5
alpha = 0.05

input_dir = Path("data/da/simulations")
output_dir = Path("results/da/clr_wilcoxon")
output_dir.mkdir(parents=True, exist_ok=True)

performance_results = []

for sim in range(1, n_simulations + 1):

    sim_dir = input_dir / f"simulation_{sim:03d}"

    # Load counts, metadata, and ground truth
    counts = pd.read_csv(
        sim_dir / "da_counts.csv", index_col=0
    )
    metadata = pd.read_csv(
        sim_dir / "da_metadata.csv"
    )
    truth = pd.read_csv(
        sim_dir / "da_ground_truth.csv"
    )

    # Match metadata to count table
    metadata = metadata.set_index("SampleID")
    metadata = metadata.loc[counts.columns]

    # CLR transformation
    log_counts = np.log(counts + pseudocount)
    clr_counts = log_counts.subtract(
        log_counts.mean(axis=0), axis=1
    )

    # Identify groups
    control_samples = metadata.index[
        metadata["Group"] == "Control"
    ]
    case_samples = metadata.index[
        metadata["Group"] == "Case"
    ]

    results = []

    # Wilcoxon rank-sum test for each taxon
    for taxon in clr_counts.index:

        control_values = clr_counts.loc[
            taxon, control_samples
        ]
        case_values = clr_counts.loc[
            taxon, case_samples
        ]

        statistic, p_value = ranksums(
            case_values, control_values
        )

        effect = (
            case_values.median()
            - control_values.median()
        )

        results.append({
            "Taxon": taxon,
            "Statistic": statistic,
            "Effect": effect,
            "PValue": p_value
        })

    results = pd.DataFrame(results)

    # Benjamini-Hochberg correction
    results["Adjusted_PValue"] = multipletests(
        results["PValue"],
        method="fdr_bh"
    )[1]

    results["Significant"] = (
        results["Adjusted_PValue"] < alpha
    )

    # Add ground truth
    results = results.merge(
        truth, on="Taxon", how="left",
        validate="one_to_one"
    )

    if results["True_DA"].isna().any():
        raise ValueError(
            f"Missing ground truth in simulation {sim}"
        )

    # Confusion matrix
    tp = int((
        results["Significant"]
        & (results["True_DA"] == 1)
    ).sum())

    fp = int((
        results["Significant"]
        & (results["True_DA"] == 0)
    ).sum())

    fn = int((
        ~results["Significant"]
        & (results["True_DA"] == 1)
    ).sum())

    tn = int((
        ~results["Significant"]
        & (results["True_DA"] == 0)
    ).sum())

    # Performance metrics
    fdp = fp / (tp + fp) if (tp + fp) > 0 else 0
    power = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    sensitivity = power
    specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan

    performance_results.append({
        "Simulation": sim,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "Discoveries": tp + fp,
        "FDP": fdp,
        "Power": power,
        "Sensitivity": sensitivity,
        "Specificity": specificity
    })

    # Save individual simulation results
    results = results.sort_values("Adjusted_PValue")
    results.to_csv(
        output_dir / f"simulation_{sim:03d}_results.csv",
        index=False
    )

    print(f"Simulation {sim}/100 analyzed.")

# Save performance across simulations
performance_df = pd.DataFrame(performance_results)

performance_df.to_csv(
    output_dir / "performance_by_simulation.csv",
    index=False
)

# Summary statistics
metrics = [
    "TP", "FP", "FN", "TN", "Discoveries",
    "FDP", "Power", "Sensitivity", "Specificity"
]

summary = performance_df[metrics].agg(
    ["mean", "std", "min", "max"]
).T

summary.to_csv(
    output_dir / "performance_summary.csv"
)

print("\nAll 100 CLR + Wilcoxon analyses completed!")
print("\nPerformance summary:")
print(summary.round(4).to_string())
print("\nResults saved to:", output_dir)
