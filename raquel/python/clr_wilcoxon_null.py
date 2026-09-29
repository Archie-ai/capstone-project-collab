import numpy as np
import pandas as pd
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import os

# Load null count data and sample metadata
counts = pd.read_csv("data/null/null_counts.csv", index_col=0)
metadata = pd.read_csv("data/null/null_metadata.csv")

# Match metadata order to the samples in the count table
metadata = metadata.set_index("SampleID")
metadata = metadata.loc[counts.columns]

# Add a small pseudocount to avoid log(0)
counts_pseudo = counts + 0.5

# CLR transformation
# For each sample, subtract the mean log abundance
log_counts = np.log(counts_pseudo)
clr_counts = log_counts.subtract(log_counts.mean(axis=0), axis=1)

# Identify Control and Case samples
control_samples = metadata.index[metadata["Group"] == "Control"]
case_samples = metadata.index[metadata["Group"] == "Case"]

results = []

# Compare Control and Case for each taxon
for taxon in clr_counts.index:

    control_values = clr_counts.loc[taxon, control_samples]
    case_values = clr_counts.loc[taxon, case_samples]

    # Wilcoxon rank-sum test
    statistic, p_value = ranksums(case_values, control_values)

    # Difference in median CLR abundance
    effect = case_values.median() - control_values.median()

    results.append({
        "Taxon": taxon,
        "Statistic": statistic,
        "Effect": effect,
        "PValue": p_value
    })

# Convert results to a table
results = pd.DataFrame(results)

# Adjust p-values using Benjamini-Hochberg FDR correction
results["Adjusted_PValue"] = multipletests(
    results["PValue"],
    method="fdr_bh"
)[1]

# Significant if adjusted p-value is below 0.05
results["Significant"] = results["Adjusted_PValue"] < 0.05

# Sort by adjusted p-value
results = results.sort_values("Adjusted_PValue")

# Create output folder and save results
os.makedirs("results/null", exist_ok=True)

results.to_csv(
    "results/null/clr_wilcoxon_results.csv",
    index=False
)

# Print a short summary
print("CLR/Wilcoxon analysis complete.")
print("Taxa tested:", len(results))
print("Significant taxa:", results["Significant"].sum())
print("Smallest adjusted p-value:", round(results["Adjusted_PValue"].min(), 4))
print("Results saved to results/null/clr_wilcoxon_results.csv")
