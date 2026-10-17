from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("erinmaley/data/null_benchmark_12")
BASE.mkdir(parents=True, exist_ok=True)

N_SIMULATIONS = 5
N_TAXA = 100
N_SAMPLES = 50

taxa = [f"Taxon_{i:03d}" for i in range(1, 101)]
samples = [f"Sample_{i:03d}" for i in range(1, 51)]
groups = ["Control"] * 25 + ["Case"] * 25

for sim in range(1, N_SIMULATIONS + 1):
    rng = np.random.default_rng(3000 + sim)

    baseline = rng.dirichlet(np.repeat(0.5, N_TAXA))
    counts = np.zeros((N_TAXA, N_SAMPLES), dtype=int)

    for j in range(N_SAMPLES):
        sample_probs = rng.dirichlet(
            np.maximum(baseline * 200, 1e-8)
        )
        library_size = max(1, rng.poisson(10000))
        counts[:, j] = rng.multinomial(
            library_size, sample_probs
        )

    folder = BASE / f"simulation_{sim:03d}"
    folder.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(
        counts, index=taxa, columns=samples
    ).rename_axis("Taxon").to_csv(folder / "counts.csv")

    pd.DataFrame({
        "SampleID": samples,
        "Group": groups
    }).to_csv(folder / "metadata.csv", index=False)

    pd.DataFrame({
        "Taxon": taxa,
        "Planted_Target": [False] * 100,
        "True_DA": [False] * 100,
        "True_Direction": ["No_change"] * 100
    }).to_csv(folder / "ground_truth.csv", index=False)

    print(f"Null simulation {sim}/5 generated")

print("NULL BENCHMARK DATASETS COMPLETE")
