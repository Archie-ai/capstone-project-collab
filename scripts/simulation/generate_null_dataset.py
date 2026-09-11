from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# CONFIGURATION
# ============================================================

RANDOM_SEED = 42

N_SAMPLES = 50
N_TAXA = 100

N_CONTROL = 25
N_CASE = 25

MEAN_LIBRARY_SIZE = 10000

# Dirichlet concentration controls how uneven the underlying
# microbial community is.
DIRICHLET_CONCENTRATION = 0.5


# ============================================================
# OUTPUT PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "simulation_null"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COUNT_FILE = OUTPUT_DIR / "null_counts.csv"
METADATA_FILE = OUTPUT_DIR / "null_metadata.csv"
TRUTH_FILE = OUTPUT_DIR / "null_ground_truth.csv"
SUMMARY_FILE = OUTPUT_DIR / "null_simulation_summary.txt"


# ============================================================
# RANDOM NUMBER GENERATOR
# ============================================================

rng = np.random.default_rng(RANDOM_SEED)


# ============================================================
# SAMPLE AND TAXON IDENTIFIERS
# ============================================================

sample_ids = [f"Sample_{i:03d}" for i in range(1, N_SAMPLES + 1)]
taxon_ids = [f"Taxon_{i:03d}" for i in range(1, N_TAXA + 1)]

groups = np.array(
    ["Control"] * N_CONTROL
    + ["Case"] * N_CASE
)


# ============================================================
# SIMULATE BASELINE MICROBIAL COMMUNITY
# ============================================================

# All samples in both groups are generated from the SAME
# underlying abundance distribution.
#
# Therefore:
# TRUE NUMBER OF DIFFERENTIALLY ABUNDANT TAXA = 0

alpha = np.repeat(DIRICHLET_CONCENTRATION, N_TAXA)

baseline_proportions = rng.dirichlet(alpha)


# ============================================================
# SIMULATE SAMPLE LIBRARY SIZES
# ============================================================

library_sizes = rng.poisson(
    lam=MEAN_LIBRARY_SIZE,
    size=N_SAMPLES
)

# Prevent extremely unlikely zero-depth samples.
library_sizes = np.maximum(library_sizes, 1)


# ============================================================
# GENERATE COUNT TABLE
# ============================================================

counts = np.zeros(
    (N_TAXA, N_SAMPLES),
    dtype=int
)

for sample_index in range(N_SAMPLES):

    counts[:, sample_index] = rng.multinomial(
        n=library_sizes[sample_index],
        pvals=baseline_proportions
    )


# ============================================================
# CREATE COUNT DATAFRAME
# ============================================================

count_df = pd.DataFrame(
    counts,
    index=taxon_ids,
    columns=sample_ids
)

count_df.index.name = "Taxon"


# ============================================================
# CREATE SAMPLE METADATA
# ============================================================

metadata_df = pd.DataFrame(
    {
        "SampleID": sample_ids,
        "Group": groups,
        "LibrarySize": library_sizes
    }
)


# ============================================================
# CREATE GROUND-TRUTH TABLE
# ============================================================

truth_df = pd.DataFrame(
    {
        "Taxon": taxon_ids,
        "True_DA": 0,
        "True_Direction": "None",
        "True_Effect_Size": 0.0
    }
)


# ============================================================
# BASIC QUALITY-CONTROL METRICS
# ============================================================

observed_library_sizes = count_df.sum(axis=0)

zero_fraction = (count_df == 0).sum().sum() / count_df.size

mean_taxon_count = count_df.to_numpy().mean()

median_taxon_count = np.median(count_df.to_numpy())


# ============================================================
# SAVE FILES
# ============================================================

count_df.to_csv(COUNT_FILE)
metadata_df.to_csv(METADATA_FILE, index=False)
truth_df.to_csv(TRUTH_FILE, index=False)


# ============================================================
# WRITE SIMULATION SUMMARY
# ============================================================

summary_text = f"""
NULL MICROBIOME SIMULATION SUMMARY
==================================

Random seed:
{RANDOM_SEED}

Samples:
{N_SAMPLES}

Control samples:
{N_CONTROL}

Case samples:
{N_CASE}

Taxa:
{N_TAXA}

True differentially abundant taxa:
0

Mean target library size:
{MEAN_LIBRARY_SIZE}

Observed mean library size:
{observed_library_sizes.mean():.2f}

Observed minimum library size:
{observed_library_sizes.min()}

Observed maximum library size:
{observed_library_sizes.max()}

Dirichlet concentration:
{DIRICHLET_CONCENTRATION}

Fraction of zero cells:
{zero_fraction:.4f}

Mean count per taxon-sample cell:
{mean_taxon_count:.2f}

Median count per taxon-sample cell:
{median_taxon_count:.2f}

Interpretation:
Case and control samples were generated from the same underlying
microbial abundance distribution. Therefore, no taxon is truly
differentially abundant in this dataset.
"""

SUMMARY_FILE.write_text(summary_text.strip() + "\n")


# ============================================================
# TERMINAL OUTPUT
# ============================================================

print("Null microbiome simulation completed successfully.")
print()
print(f"Count table:       {COUNT_FILE}")
print(f"Metadata:          {METADATA_FILE}")
print(f"Ground truth:      {TRUTH_FILE}")
print(f"Simulation summary:{SUMMARY_FILE}")
print()
print(f"Count table shape: {count_df.shape}")
print(f"Zero fraction:     {zero_fraction:.4f}")
print(
    "Mean library size:"
    f" {observed_library_sizes.mean():.2f}"
)
