#!/usr/bin/env python3

# ============================================================
# CDI SCHUBERT REAL-DATA FORMATTING
# ============================================================
#
# This script converts the original MicrobiomeHD cdi_schubert
# files into the count-table and metadata format used by the
# consensus pipeline.
#
# Original source files are not changed.
#
# This script was developed with assistance from OpenAI's ChatGPT
#

from pathlib import Path

import pandas as pd


# ============================================================
# FILE LOCATIONS
# ============================================================

metadata_file = Path(
    "erinmaley/data/real_test_input/raw_source/extracted/"
    "cdi_schubert_results/cdi_schubert.metadata.txt"
)

otu_file = Path(
    "erinmaley/data/real_test_input/raw_source/extracted/"
    "cdi_schubert_results/RDP/"
    "cdi_schubert.otu_table.100.denovo.rdp_assigned"
)

output_dir = Path(
    "erinmaley/data/real_test_input/formatted"
)

output_dir.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# LOAD SOURCE DATA
# ============================================================

# The original metadata contains characters that are not UTF-8,
# so latin1 is used only when reading the source file.

metadata = pd.read_csv(
    metadata_file,
    sep="\t",
    encoding="latin1"
)

otu = pd.read_csv(
    otu_file,
    sep="\t",
    encoding="latin1"
)


# ============================================================
# SELECT THE BINARY COMPARISON
# ============================================================

# This pipeline currently compares two groups:
#
# H   -> Control
# CDI -> Case
#
# nonCDI samples are not included in this run.

metadata = metadata[
    metadata["DiseaseState"].isin(
        ["H", "CDI"]
    )
].copy()

metadata["sample_id"] = (
    metadata["sample_id"]
    .astype(str)
    .str.strip()
)


# ============================================================
# KEEP ONLY SAMPLES PRESENT IN BOTH FILES
# ============================================================

otu_sample_ids = set(
    otu.columns[1:].astype(str)
)

metadata = metadata[
    metadata["sample_id"].isin(
        otu_sample_ids
    )
].copy()


# ============================================================
# CREATE PIPELINE METADATA
# ============================================================

group_map = {
    "H": "Control",
    "CDI": "Case"
}

formatted_metadata = pd.DataFrame({
    "SampleID": metadata["sample_id"],
    "Group": metadata["DiseaseState"].map(group_map)
})

# Remove accidental duplicates if the source metadata ever
# contains the same sample more than once.
formatted_metadata = (
    formatted_metadata
    .drop_duplicates(subset="SampleID")
    .reset_index(drop=True)
)


# ============================================================
# CREATE PIPELINE COUNT TABLE
# ============================================================

sample_ids = formatted_metadata[
    "SampleID"
].tolist()

# The first OTU-table column contains the taxonomic/OTU label.
taxon_column = otu.columns[0]

formatted_counts = otu[
    [taxon_column] + sample_ids
].copy()

formatted_counts = formatted_counts.rename(
    columns={
        taxon_column: "Taxon"
    }
)


# ============================================================
# VERIFY COUNTS
# ============================================================

count_values = formatted_counts.iloc[:, 1:]

# All abundance values should be numeric.
count_values = count_values.apply(
    pd.to_numeric,
    errors="raise"
)

# The pipeline expects non-negative integer counts.
if (count_values < 0).any().any():
    raise ValueError(
        "Negative abundance values found."
    )

if not (
    count_values
    == count_values.round()
).all().all():
    raise ValueError(
        "Non-integer abundance values found."
    )

formatted_counts.iloc[:, 1:] = (
    count_values.astype("int64")
)


# ============================================================
# FINAL SAMPLE CHECK
# ============================================================

count_samples = set(
    formatted_counts.columns[1:]
)

metadata_samples = set(
    formatted_metadata["SampleID"]
)

if count_samples != metadata_samples:
    raise ValueError(
        "Formatted count and metadata sample IDs do not match."
    )


# ============================================================
# SAVE FORMATTED INPUTS
# ============================================================

counts_output = (
    output_dir
    / "cdi_schubert_counts.csv"
)

metadata_output = (
    output_dir
    / "cdi_schubert_metadata.csv"
)

formatted_counts.to_csv(
    counts_output,
    index=False
)

formatted_metadata.to_csv(
    metadata_output,
    index=False
)


# ============================================================
# SUMMARY
# ============================================================

print()
print("CDI SCHUBERT FORMATTING COMPLETED")
print("=================================")

print(
    f"Samples retained: "
    f"{len(formatted_metadata)}"
)

print(
    f"Control: "
    f"{(formatted_metadata['Group'] == 'Control').sum()}"
)

print(
    f"Case: "
    f"{(formatted_metadata['Group'] == 'Case').sum()}"
)

print(
    f"Taxa/OTUs: "
    f"{formatted_counts.shape[0]}"
)

print()
print("Formatted count table:")
print(f"  {counts_output}")

print()
print("Formatted metadata:")
print(f"  {metadata_output}")

print()
