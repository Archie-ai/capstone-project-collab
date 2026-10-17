from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("erinmaley/data/null_benchmark_12")
df = pd.read_csv(BASE / "null_200_method_summary.csv")

methods = ["ALDEx2", "ANCOM_BC", "CLR_Wilcoxon", "DESeq2"]
labels = ["ALDEx2", "ANCOM-BC", "CLR/Wilcoxon", "DESeq2"]
colors = ["#293B69", "#7189B6", "#B28AA8", "#D49A76"]

df = df.set_index("Method").loc[methods]
x = np.arange(4)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False
})

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

metrics = [
    ("MeanTypeIError", "A. Empirical Per-Taxon Type I Error", 100),
    ("ProbabilityAnyFalseDiscovery",
     "B. Simulations With ≥1 False Discovery", 100)
]

for ax, (column, title, multiplier) in zip(axes, metrics):
    values = df[column].to_numpy() * multiplier

    bars = ax.bar(x, values, color=colors, width=0.6)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width()/2,
            value + max(values)*0.04,
            f"{value:.3f}%" if column == "MeanTypeIError"
            else f"{value:.1f}%",
            ha="center", fontsize=9
        )

    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel("Percentage (%)")
    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_ylim(0, max(values)*1.25)
    ax.grid(axis="y", alpha=0.15)
    ax.set_axisbelow(True)

fig.suptitle(
    "Empirical False-Positive Behavior Across 200 Null Simulations",
    fontsize=13, weight="bold"
)

fig.tight_layout()

output = BASE / "figure4_null_200_method_performance.png"
fig.savefig(output, dpi=300, bbox_inches="tight")
plt.close(fig)

print("200-SIMULATION FIGURE COMPLETE:", output)
