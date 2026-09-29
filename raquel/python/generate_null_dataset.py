import numpy as np
import pandas as pd
import os

# Simulation settings
np.random.seed(42)

n_samples = 50
n_taxa = 100
mean_library_size = 10000

# Create sample and taxon names
sample_ids = [f"Sample_{i+1}" for i in range(n_samples)]
taxa_ids = [f"Taxon_{i+1}" for i in range(n_taxa)]

# 25 control samples and 25 case samples
groups = ["Control"] * 25 + ["Case"] * 25

# Create one baseline microbial composition
# Both groups use the same composition, so there are no true DA taxa
alpha = np.repeat(0.5, n_taxa)
baseline = np.random.dirichlet(alpha)

# Simulate sequencing depth for each sample
library_sizes = np.random.poisson(mean_library_size, n_samples)

# Generate count data
counts = []

for library_size in library_sizes:
    sample_counts = np.random.multinomial(library_size, baseline)
    counts.append(sample_counts)

# Convert to taxa x samples format
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

# Create ground-truth table
# Every taxon is non-DA in this null simulation
truth_df = pd.DataFrame({
    "Taxon": taxa_ids,
    "True_DA": 0
})

# Create output folder
os.makedirs("data/null", exist_ok=True)

# Save files
counts_df.to_csv("data/null/null_counts.csv")
metadata_df.to_csv("data/null/null_metadata.csv", index=False)
truth_df.to_csv("data/null/null_ground_truth.csv", index=False)

# Basic QC
zero_fraction = (counts == 0).mean()

print("Null dataset generated successfully.")
print("Samples:", n_samples)
print("Taxa:", n_taxa)
print("True DA taxa: 0")
print("Mean library size:", round(library_sizes.mean(), 2))
print("Zero fraction:", round(zero_fraction, 4))
