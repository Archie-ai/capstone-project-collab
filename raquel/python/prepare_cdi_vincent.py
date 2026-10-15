# BIOT 670I - Microbiome Abundance Consensus Machine
# Prepare the cdi_vincent real microbiome dataset

from pathlib import Path
import pandas as pd

# File locations
RAW = Path("data/raw/cdi_vincent")
OUTPUT = Path("data/processed/cdi_vincent")
OUTPUT.mkdir(parents=True, exist_ok=True)

COUNTS_FILE = RAW / (
    "cdi_vincent_v3v5.otu_table.100.denovo.rdp_assigned"
)
METADATA_FILE = RAW / "cdi_vincent_v3v5.metadata.txt"

# Minimum prevalence: 5 of 50 samples (10%)
MIN_PREVALENCE = 0.10

# Load original files
counts = pd.read_csv(
    COUNTS_FILE,
    sep="\t",
    index_col=0
)

metadata = pd.read_csv(
    METADATA_FILE,
    sep="\t",
    index_col=0
)

# Validate sample IDs
if counts.columns.duplicated().any():
    raise ValueError("Duplicate sample IDs in OTU table")

if metadata.index.duplicated().any():
    raise ValueError("Duplicate sample IDs in metadata")

if set(counts.columns) != set(metadata.index):
    raise ValueError("OTU table and metadata samples do not match")

# Ensure samples are in the same order
metadata = metadata.loc[counts.columns]

# Define Case and Control groups
group_mapping = {
    "CDI": "Case",
    "H": "Control"
}

if not metadata["DiseaseState"].isin(group_mapping).all():
    raise ValueError("Unexpected DiseaseState category")

metadata["Group"] = metadata["DiseaseState"].map(group_mapping)

# Validate counts
counts = counts.apply(pd.to_numeric, errors="raise")

if counts.isna().any().any():
    raise ValueError("Missing count values detected")

if (counts < 0).any().any():
    raise ValueError("Negative counts detected")

if ((counts % 1) != 0).any().any():
    raise ValueError("Non-integer counts detected")

counts = counts.astype("int64")

# Calculate prevalence using all 50 samples
prevalence_count = (counts > 0).sum(axis=1)

minimum_samples = max(
    1,
    int(MIN_PREVALENCE * counts.shape[1] + 0.999999)
)

# Filter rare OTUs
keep = prevalence_count >= minimum_samples
filtered_counts = counts.loc[keep].copy()

# Preserve taxonomy labels and generate unique OTU IDs
otu_ids = [
    f"OTU_{i:04d}"
    for i in range(1, len(filtered_counts) + 1)
]

taxonomy_mapping = pd.DataFrame({
    "OTU_ID": otu_ids,
    "Taxonomy": filtered_counts.index,
    "Prevalence_Count": prevalence_count.loc[keep].values,
    "Prevalence_Fraction": (
        prevalence_count.loc[keep].values / counts.shape[1]
    )
})

filtered_counts.index = otu_ids
filtered_counts.index.name = "Taxon"

# Standardized metadata
processed_metadata = pd.DataFrame({
    "SampleID": metadata.index,
    "Group": metadata["Group"].values,
    "DiseaseState": metadata["DiseaseState"].values
})

# Save processed data
filtered_counts.to_csv(
    OUTPUT / "cdi_counts.csv"
)

processed_metadata.to_csv(
    OUTPUT / "cdi_metadata.csv",
    index=False
)

taxonomy_mapping.to_csv(
    OUTPUT / "cdi_taxonomy_mapping.csv",
    index=False
)

# Save preprocessing summary
summary = pd.DataFrame([{
    "Original_OTUs": len(counts),
    "Retained_OTUs": len(filtered_counts),
    "Removed_OTUs": len(counts) - len(filtered_counts),
    "Samples": counts.shape[1],
    "Case_Samples": int((metadata["Group"] == "Case").sum()),
    "Control_Samples": int((metadata["Group"] == "Control").sum()),
    "Minimum_Prevalence": MIN_PREVALENCE,
    "Minimum_Samples": minimum_samples,
    "Original_Zero_Fraction": (counts == 0).to_numpy().mean(),
    "Filtered_Zero_Fraction": (
        filtered_counts == 0
    ).to_numpy().mean()
}])

summary.to_csv(
    OUTPUT / "preprocessing_summary.csv",
    index=False
)

print("cdi_vincent preprocessing completed!")
print(f"Original OTUs: {len(counts)}")
print(f"Retained OTUs: {len(filtered_counts)}")
print(f"Removed OTUs: {len(counts) - len(filtered_counts)}")
print(f"Samples: {counts.shape[1]}")
print(f"Case samples: {(metadata['Group'] == 'Case').sum()}")
print(f"Control samples: {(metadata['Group'] == 'Control').sum()}")
print(f"Minimum prevalence: {minimum_samples} samples")
print(f"Original zero fraction: {(counts == 0).to_numpy().mean():.4f}")
print(
    "Filtered zero fraction: "
    f"{(filtered_counts == 0).to_numpy().mean():.4f}"
)
print(f"Files saved to: {OUTPUT}")
