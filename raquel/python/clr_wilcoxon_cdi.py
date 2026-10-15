# BIOT 670I - Microbiome Abundance Consensus Machine
# Real-data differential abundance analysis: CLR + Wilcoxon
# Dataset: cdi_vincent (CDI vs Healthy)

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

# File locations
DATA = Path("data/processed/cdi_vincent")
OUTPUT = Path("results/real/cdi_vincent/clr_wilcoxon")
OUTPUT.mkdir(parents=True, exist_ok=True)

# Analysis settings
PSEUDOCOUNT = 0.5
ALPHA = 0.05

# Load processed dataset
counts = pd.read_csv(
    DATA / "cdi_counts_filtered.csv",
    index_col=0
)

metadata = pd.read_csv(
    DATA / "cdi_metadata.csv"
)

taxonomy = pd.read_csv(
    DATA / "cdi_taxonomy_mapping.csv"
)

# Validate and align samples
if counts.index.duplicated().any():
    raise ValueError("Duplicate OTU IDs")

if counts.columns.duplicated().any():
    raise ValueError("Duplicate sample IDs")

if metadata["SampleID"].duplicated().any():
    raise ValueError("Duplicate metadata sample IDs")

if set(counts.columns) != set(metadata["SampleID"]):
    raise ValueError("Count and metadata sample IDs do not match")

metadata = metadata.set_index("SampleID").loc[counts.columns]

if not metadata["Group"].isin(["Case", "Control"]).all():
    raise ValueError("Unexpected sample group")

# Samples x OTUs
sample_counts = counts.T.astype(float)

# CLR transformation
log_counts = np.log(sample_counts + PSEUDOCOUNT)

clr = log_counts.sub(
    log_counts.mean(axis=1),
    axis=0
)

# Separate CDI and healthy groups
case_ids = metadata.index[metadata["Group"] == "Case"]
control_ids = metadata.index[metadata["Group"] == "Control"]

results = []

for otu in clr.columns:

    case_values = clr.loc[case_ids, otu]
    control_values = clr.loc[control_ids, otu]

    statistic, p_value = ranksums(
        case_values,
        control_values
    )

    # Positive effect means higher CLR abundance in CDI
    effect = case_values.mean() - control_values.mean()

    results.append({
        "Taxon": otu,
        "Statistic": statistic,
        "Effect": effect,
        "PValue": p_value
    })

results = pd.DataFrame(results)

# Benjamini-Hochberg adjustment
results["Adjusted_PValue"] = multipletests(
    results["PValue"],
    method="fdr_bh"
)[1]

results["Significant"] = (
    results["Adjusted_PValue"] < ALPHA
)

# Attach taxonomy information
results = results.merge(
    taxonomy[["OTU_ID", "Taxonomy"]],
    left_on="Taxon",
    right_on="OTU_ID",
    how="left",
    validate="one_to_one"
).drop(columns="OTU_ID")

# Sort by adjusted p-value
results = results.sort_values(
    "Adjusted_PValue"
).reset_index(drop=True)

# Save results
results.to_csv(
    OUTPUT / "clr_wilcoxon_results.csv",
    index=False
)

# Save CLR-transformed abundance table
clr.T.to_csv(
    OUTPUT / "clr_transformed_counts.csv",
    index_label="Taxon"
)

# Summary
summary = pd.DataFrame([{
    "Method": "CLR_Wilcoxon",
    "Dataset": "cdi_vincent",
    "Samples": len(metadata),
    "Case_Samples": len(case_ids),
    "Control_Samples": len(control_ids),
    "OTUs_Tested": len(results),
    "Significant_OTUs": int(results["Significant"].sum()),
    "Alpha": ALPHA,
    "Pseudocount": PSEUDOCOUNT
}])

summary.to_csv(
    OUTPUT / "analysis_summary.csv",
    index=False
)

print("CLR + Wilcoxon analysis completed!")
print(f"Samples: {len(metadata)}")
print(f"CDI samples: {len(case_ids)}")
print(f"Healthy samples: {len(control_ids)}")
print(f"OTUs tested: {len(results)}")
print(f"Significant OTUs: {results['Significant'].sum()}")
print(f"Minimum adjusted p-value: {results['Adjusted_PValue'].min():.6g}")
print(f"Results saved to: {OUTPUT}")
