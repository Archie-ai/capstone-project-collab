from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("erinmaley/data/benchmark_11")
method = pd.read_csv(BASE / "scored_summary.csv")
consensus = pd.read_csv(BASE / "consensus_scored_summary.csv")

names = ["ALDEx2", "ANCOM_BC", "CLR_Wilcoxon", "DESeq2"]
labels = ["ALDEx2", "ANCOM-BC", "CLR/Wilcoxon", "DESeq2"]
colors = ["#293B69", "#7189B6", "#B28AA8", "#D49A76"]

method = method.set_index("Method").loc[names]
consensus = consensus.sort_values("Threshold")

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False
})

def panel(ax, labels, values, errors, colors, title):
    x = np.arange(len(labels))
    ax.bar(x, values * 100, color=colors, width=0.62)
    ax.errorbar(
        x, values * 100, yerr=errors * 100,
        fmt="none", ecolor="#333333",
        capsize=4, linewidth=1.1
    )
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel("Percentage (%)")
    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_ylim(0, max(5, np.max((values + errors) * 100) * 1.2))
    ax.grid(axis="y", alpha=0.15)
    ax.set_axisbelow(True)

# Figure 1: Individual method performance
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
panel(
    axes[0], labels,
    method["MeanFDP"].to_numpy(),
    method["SDFDP"].to_numpy(),
    colors, "A. Mean planted-target FDP"
)
panel(
    axes[1], labels,
    method["MeanSensitivity"].to_numpy(),
    method["SDSensitivity"].to_numpy(),
    colors, "B. Mean planted-target sensitivity"
)
fig.suptitle(
    "Differential Abundance Method Performance",
    fontsize=13, weight="bold"
)
fig.tight_layout()
fig.savefig(BASE / "figure1_method_performance.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)

# Figure 2: Consensus performance
votes = [f"{int(n)} of 4" for n in consensus["Threshold"]]
vote_colors = ["#293B69", "#B28AA8", "#D49A76"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
panel(
    axes[0], votes,
    consensus["MeanFDP"].to_numpy(),
    consensus["SDFDP"].to_numpy(),
    vote_colors, "A. Mean planted-target FDP"
)
panel(
    axes[1], votes,
    consensus["MeanSensitivity"].to_numpy(),
    consensus["SDSensitivity"].to_numpy(),
    vote_colors, "B. Mean planted-target sensitivity"
)

axes[0].text(
    0.5, 0.85,
    "No non-target calls at\n3-of-4 or 4-of-4",
    transform=axes[0].transAxes,
    ha="center", fontsize=9,
    color="#293B69"
)

fig.suptitle(
    "Directional Consensus Voting Performance",
    fontsize=13, weight="bold"
)
fig.tight_layout()
fig.savefig(BASE / "figure2_consensus_performance.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)

print("TWO BENCHMARK 11 FIGURES COMPLETED")
print("Results: exploratory, five independent simulations")
print("Error bars: standard deviation")
