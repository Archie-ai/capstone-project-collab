from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd

BASE = Path("erinmaley/data/null_benchmark_12")
ROOT = Path("erinmaley")
TOTAL = 200

taxa = [f"Taxon_{i:03d}" for i in range(1, 101)]
samples = [f"Sample_{i:03d}" for i in range(1, 51)]

for i in range(1, TOTAL + 1):
    folder = BASE / f"simulation_{i:03d}"
    done = folder / "completed.flag"

    if done.exists():
        print(f"Simulation {i}/200 already complete", flush=True)
        continue

    folder.mkdir(parents=True, exist_ok=True)

    # Generate missing simulation inputs.
    if not (folder / "counts.csv").exists():
        rng = np.random.default_rng(3000 + i)
        baseline = rng.dirichlet(np.repeat(0.5, 100))
        counts = np.zeros((100, 50), dtype=int)

        for j in range(50):
            probabilities = rng.dirichlet(
                np.maximum(baseline * 200, 1e-8)
            )
            depth = max(1, rng.poisson(10000))
            counts[:, j] = rng.multinomial(depth, probabilities)

        pd.DataFrame(
            counts, index=taxa, columns=samples
        ).rename_axis("Taxon").to_csv(folder / "counts.csv")

        pd.DataFrame({
            "SampleID": samples,
            "Group": ["Control"] * 25 + ["Case"] * 25
        }).to_csv(folder / "metadata.csv", index=False)

        pd.DataFrame({
            "Taxon": taxa,
            "Planted_Target": [False] * 100,
            "True_DA": [False] * 100,
            "True_Direction": ["No_change"] * 100
        }).to_csv(folder / "ground_truth.csv", index=False)

    qc = folder / "qc"
    out = folder / "method_outputs"

    commands = [
        [
            "python3", str(ROOT / "scripts/quality_control.py"),
            "--counts", str(folder / "counts.csv"),
            "--metadata", str(folder / "metadata.csv"),
            "--output-dir", str(qc),
            "--dataset-label", "signal"
        ]
    ]

    for method, script in [
        ("aldex2", "aldex2.R"),
        ("ancombc", "ancombc.R"),
        ("deseq2", "deseq2.R")
    ]:
        commands.append([
            "Rscript", str(ROOT / "scripts/methods_3" / script),
            str(qc / "qc_signal_counts.csv"),
            str(qc / "qc_signal_metadata.csv"),
            str(out / method),
            "signal"
        ])

    commands.append([
        "python3", str(ROOT / "scripts/methods_3/clr_wilcoxon.py"),
        "--counts", str(qc / "qc_signal_counts.csv"),
        "--metadata", str(qc / "qc_signal_metadata.csv"),
        "--output-dir", str(out / "clr_wilcoxon"),
        "--dataset-label", "signal"
    ])

    start = time.time()

    print(f"Starting simulation {i}/200", flush=True)

    try:
        with open(folder / "run.log", "a") as log:
            for command in commands:
                subprocess.run(
                    command,
                    check=True,
                    stdout=log,
                    stderr=subprocess.STDOUT
                )

        done.write_text("complete\n")

        print(
            f"Completed simulation {i}/200 "
            f"in {(time.time() - start):.1f} seconds",
            flush=True
        )

    except subprocess.CalledProcessError as exc:
        print(f"ERROR simulation {i}: {exc}", flush=True)
        raise

print("ALL 200 NULL SIMULATIONS COMPLETE")
