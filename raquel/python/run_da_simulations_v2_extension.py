import numpy as np
import pandas as pd
from pathlib import Path

# Simulation settings
n_simulations = 100
n_samples = 50
n_taxa = 100
mean_library_size = 10000

# Differential abundance settings
n_da_taxa = 10
effect_size = 2.0

# Controls sample-to-sample biological variability
# Smaller values produce greater variability
concentration = 200

# Output directory (preserves Version 1)
output_dir = Path("data/da/simulations_v2")
output_dir.mkdir(parents=True, exist_ok=True)

sample_ids = [f"Sample_{i+1}" for i in range(n_samples)]
taxa_ids = [f"Taxon_{i+1}" for i in range(n_taxa)]

groups = ["Control"] * 25 + ["Case"] * 25

qc_results = []

for sim in range(101, 201):

    seed = 41 + sim
    rng = np.random.default_rng(seed)

    # Generate baseline microbial composition
    alpha = np.repeat(0.5, n_taxa)
    baseline = rng.dirichlet(alpha)

    # Define group-level compositions
    control_composition = baseline.copy()
    case_composition = baseline.copy()

    # Introduce DA in the first 10 taxa
    case_composition[:n_da_taxa] *= effect_size
    case_composition /= case_composition.sum()

    # Generate sequencing depths
    library_sizes = rng.poisson(
        mean_library_size, n_samples
    )

    counts = []

    for i, library_size in enumerate(library_sizes):

        # Select group-level composition
        if groups[i] == "Control":
            group_composition = control_composition
        else:
            group_composition = case_composition

        # Generate sample-specific composition
        sample_composition = rng.dirichlet(
            group_composition * concentration
        )

        # Generate observed sequencing counts
        sample_counts = rng.multinomial(
            library_size, sample_composition
        )

        counts.append(sample_counts)

    # Convert to taxa x samples
    counts = np.array(counts).T

    counts_df = pd.DataFrame(
        counts,
        index=taxa_ids,
        columns=sample_ids
    )
    counts_df.index.name = "Taxon"

    metadata_df = pd.DataFrame({
        "SampleID": sample_ids,
        "Group": groups,
        "LibrarySize": library_sizes
    })

    truth_df = pd.DataFrame({
        "Taxon": taxa_ids,
        "True_DA": [
            1 if i < n_da_taxa else 0
            for i in range(n_taxa)
        ]
    })

    sim_dir = output_dir / f"simulation_{sim:03d}"
    sim_dir.mkdir(parents=True, exist_ok=True)

    counts_df.to_csv(sim_dir / "da_counts.csv")

    metadata_df.to_csv(
        sim_dir / "da_metadata.csv",
        index=False
    )

    truth_df.to_csv(
        sim_dir / "da_ground_truth.csv",
        index=False
    )

    zero_fraction = (counts == 0).mean()

    qc_results.append({
        "Simulation": sim,
        "Seed": seed,
        "Samples": n_samples,
        "Taxa": n_taxa,
        "True_DA_Taxa": n_da_taxa,
        "Effect_Size": effect_size,
        "Concentration": concentration,
        "Mean_Library_Size": library_sizes.mean(),
        "Zero_Fraction": zero_fraction
    })

    print(f"Simulation {sim}/100 completed.")

qc_df = pd.DataFrame(qc_results)

qc_df.to_csv(
    output_dir / "simulation_qc_summary_101_200.csv",
    index=False
)

print("\nSimulations 101–200 completed!")
print("Results saved to:", output_dir)
