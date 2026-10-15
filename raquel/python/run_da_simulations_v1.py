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

# Output directory
output_dir = Path("data/da/simulations")
output_dir.mkdir(parents=True, exist_ok=True)

# Sample and taxon names
sample_ids = [f"Sample_{i+1}" for i in range(n_samples)]
taxa_ids = [f"Taxon_{i+1}" for i in range(n_taxa)]

groups = ["Control"] * 25 + ["Case"] * 25

# Store QC results
qc_results = []

for sim in range(1, n_simulations + 1):

    # Different seed for each simulation
    seed = 41 + sim
    rng = np.random.default_rng(seed)

    # Create baseline composition
    alpha = np.repeat(0.5, n_taxa)
    control_composition = rng.dirichlet(alpha)

    # Create Case composition
    case_composition = control_composition.copy()
    case_composition[:n_da_taxa] *= effect_size

    # Normalize proportions
    case_composition /= case_composition.sum()

    # Generate library sizes
    library_sizes = rng.poisson(
        mean_library_size, n_samples
    )

    # Generate count data
    counts = []

    for i, library_size in enumerate(library_sizes):

        if groups[i] == "Control":
            composition = control_composition
        else:
            composition = case_composition

        sample_counts = rng.multinomial(
            library_size, composition
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

    # Create metadata
    metadata_df = pd.DataFrame({
        "SampleID": sample_ids,
        "Group": groups,
        "LibrarySize": library_sizes
    })

    # Create ground truth
    truth_df = pd.DataFrame({
        "Taxon": taxa_ids,
        "True_DA": [
            1 if i < n_da_taxa else 0
            for i in range(n_taxa)
        ]
    })

    # Create folder for each simulation
    sim_dir = output_dir / f"simulation_{sim:03d}"
    sim_dir.mkdir(parents=True, exist_ok=True)

    # Save datasets
    counts_df.to_csv(sim_dir / "da_counts.csv")
    metadata_df.to_csv(
        sim_dir / "da_metadata.csv", index=False
    )
    truth_df.to_csv(
        sim_dir / "da_ground_truth.csv", index=False
    )

    # Quality control
    zero_fraction = (counts == 0).mean()

    qc_results.append({
        "Simulation": sim,
        "Seed": seed,
        "Samples": n_samples,
        "Taxa": n_taxa,
        "True_DA_Taxa": n_da_taxa,
        "Effect_Size": effect_size,
        "Mean_Library_Size": library_sizes.mean(),
        "Zero_Fraction": zero_fraction
    })

    print(f"Simulation {sim}/100 completed.")

# Save QC summary
qc_df = pd.DataFrame(qc_results)
qc_df.to_csv(
    output_dir / "simulation_qc_summary.csv",
    index=False
)

print("\nAll 100 DA simulations completed successfully!")
print("Results saved to:", output_dir)
