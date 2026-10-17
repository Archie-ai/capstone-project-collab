from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("erinmaley/data/null_benchmark_12")
df = pd.read_csv(BASE / "null_200_consensus_summary.csv")
df = df.sort_values("Threshold")

labels = [f"{int(n)} of 4" for n in df["Threshold"]]
colors = ["#293B69", "#B28AA8", "#D49A76"]
x = np.arange(len(df))

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False
})

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

metrics = [
    ("MeanTypeIError", "A. Empirical Per-Taxon Type I Error"),
    ("ProbabilityAnyFalseDiscovery",
     "B. Simulations With ≥1 False Discovery")
]

for ax, (column, title) in zip(axes, metrics):
    values = df[column].to_numpy() * 100

    bars = ax.bar(x, values, color=colors, width=0.6)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + max(values) * 0.05,
            f"{value:.3f}%" if column == "MeanTypeIError"
            else f"{value:.1f}%",
            ha="center", fontsize=10
        )

    ax.set_xticks(x, labels)
    ax.set_ylabel("Percentage (%)")
    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_ylim(0, max(values) * 1.3)
    ax.grid(axis="y", alpha=0.15)
    ax.set_axisbelow(True)

fig.suptitle(
    "Directional Consensus Performance Across 200 Null Simulations",
    fontsize=13, weight="bold"
)

fig.tight_layout()

output = BASE / "figure5_null_200_consensus_performance.png"
fig.savefig(output, dpi=300, bbox_inches="tight")
plt.close(fig)

print("CONSENSUS FIGURE COMPLETE:", output)
