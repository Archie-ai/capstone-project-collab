from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("erinmaley/data/benchmark_11")
N_SIMULATIONS = 5
N_SAMPLES = 50
N_TAXA = 100
N_TARGETS = 20
LIBRARY_SIZE = 10000

taxa = [f"Taxon_{i:03d}" for i in range(1, N_TAXA + 1)]
samples = [f"Sample_{i:03d}" for i in range(1, N_SAMPLES + 1)]
groups = np.array(["Control"] * 25 + ["Case"] * 25)

BASE.mkdir(parents=True, exist_ok=True)

for sim in range(1, N_SIMULATIONS + 1):
    rng = np.random.default_rng(2000 + sim)

    # Generate a shared baseline microbial composition.
    baseline = rng.dirichlet(np.repeat(0.5, N_TAXA))

    # Randomly select planted targets.
    targets = rng.choice(N_TAXA, N_TARGETS, replace=False)

    increased = targets[:10]
    decreased = targets[10:]

    # Apply absolute-scale fold changes.
    weights = baseline.copy()
    weights[increased] *= 2.0
    weights[decreased] *= 0.5

    case_probs = weights / weights.sum()

    # Generate independent biological samples.
    counts = np.zeros((N_TAXA, N_SAMPLES), dtype=int)

    for j in range(N_SAMPLES):
        reference = baseline if j < 25 else case_probs

        sample_probs = rng.dirichlet(
            np.maximum(reference * 200, 1e-8)
        )

        counts[:, j] = rng.multinomial(
            max(1, rng.poisson(LIBRARY_SIZE)),
            sample_probs
        )

    output = BASE / f"simulation_{sim:03d}"
    output.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(
        counts, index=taxa, columns=samples
    ).rename_axis("Taxon").to_csv(output / "counts.csv")

    pd.DataFrame({
        "SampleID": samples,
        "Group": groups
    }).to_csv(output / "metadata.csv", index=False)

    relative_fc = case_probs / baseline

    truth = pd.DataFrame({
        "Taxon": taxa,
        "Planted_Target": np.isin(np.arange(N_TAXA), targets),
        "True_DA": np.isin(np.arange(N_TAXA), targets),
        "Baseline_Proportion": baseline,
        "Case_Proportion": case_probs,
        "Relative_FoldChange": relative_fc,
        "True_Direction": np.where(
            relative_fc > 1,
            "Higher_in_Case",
            "Higher_in_Control"
        )
    })

    truth.to_csv(output / "ground_truth.csv", index=False)

    assert len(truth) == N_TAXA
    assert truth["Planted_Target"].sum() == N_TARGETS
    assert (counts >= 0).all()
    assert (counts.sum(axis=0) > 0).all()

    print(f"Simulation {sim}: 50 samples, 100 taxa, 20 planted targets")

print("BENCHMARK 11 PILOT GENERATED SUCCESSFULLY")
print("Note: non-target relative abundances may shift due to closure.")
