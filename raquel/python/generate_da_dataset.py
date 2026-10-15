import numpy as np
import pandas as pd
import os

# Simulation settings
np.random.seed(42)

n_samples = 50
n_taxa = 100
mean_library_size = 10000

# DA settings
n_da_taxa = 10
effect_size = 2.0

# Create sample and taxon names
sample_ids = [f"Sample_{i+1}" for i in range(n_samples)]
taxa_ids = [f"Taxon_{i+1}" for i in range(n_taxa)]

# Create groups
groups = ["Control"] * 25 + ["Case"] * 25

# Create baseline composition
alpha = np.repeat(0.5, n_taxa)
control_composition = np.random.dirichlet(alpha)

# Create Case composition
case_composition = control_composition.copy()

# Increase first 10 taxa
case_composition[:n_da_taxa] *= effect_size

# Normalize proportions
case_composition /= case_composition.sum()

# Simulate library sizes
library_sizes = np.random.poisson(mean_library_size, n_samples)

# Generate count data
counts = []

for i, library_size in enumerate(library_sizes):
    if groups[i] == "Control":
        composition = control_composition
    else:
        composition = case_composition

    sample_counts = np.random.multinomial(library_size, composition)
    counts.append(sample_counts)

# Convert to taxa x samples
counts = np.array(counts).T

counts_df = pd.DataFrame(
    counts,
    index=taxa_ids,
    columns=sample_ids
)

counts_df.index.name = "Taxon"

# Create metadata
metadata_df = pd.DataFrame({
    "SampleID": sample_ids,
    "Group": groups,
    "LibrarySize": library_sizes
})

# Create ground truth
truth_df = pd.DataFrame({
    "Taxon": taxa_ids,
    "True_DA": [1 if i < n_da_taxa else 0 for i in range(n_taxa)]
})

# Create output folder
os.makedirs("data/da", exist_ok=True)

# Save files
counts_df.to_csv("data/da/da_counts.csv")
metadata_df.to_csv("data/da/da_metadata.csv", index=False)
truth_df.to_csv("data/da/da_ground_truth.csv", index=False)

# Basic QC
zero_fraction = (counts == 0).mean()

print("DA dataset generated successfully.")
print("Samples:", n_samples)
print("Taxa:", n_taxa)
print("True DA taxa:", truth_df["True_DA"].sum())
print("Effect size:", effect_size)
print("Mean library size:", round(library_sizes.mean(), 2))
print("Zero fraction:", round(zero_fraction, 4))
