# BIOT 670I - Microbiome Abundance Consensus Machine
# Figure 2: Consensus Voting Performance Across 200 Simulations

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Input and output directories
BASE = Path("results/da/consensus_v2_200")
OUTPUT = Path("results/figures")
OUTPUT.mkdir(parents=True, exist_ok=True)

# Load validated consensus performance results
performance = pd.read_csv(
    BASE / "consensus_performance_by_simulation.csv"
)

# Voting rules
RULES = ["2-of-4", "3-of-4", "4-of-4"]

# Calculate mean performance for each consensus rule
summary = (
    performance.groupby("Rule")[["FDP", "Power"]]
    .mean()
    .reindex(RULES)
    .reset_index()
)

print("\nConsensus Voting Performance (200 Simulations)")
print(summary.to_string(index=False))

# Create a two-panel comparison figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

colors = ["#52796F", "#84A98C", "#D4A373"]

# Panel A: Mean false discovery proportion
bars1 = axes[0].bar(
    summary["Rule"],
    summary["FDP"],
    color=colors,
    edgecolor="black",
    linewidth=0.5
)

axes[0].set_title("A. False Discovery Proportion", fontweight="bold")
axes[0].set_ylabel("Mean FDP")
axes[0].set_ylim(0, 0.06)

# Panel B: Mean statistical power
bars2 = axes[1].bar(
    summary["Rule"],
    summary["Power"],
    color=colors,
    edgecolor="black",
    linewidth=0.5
)

axes[1].set_title("B. Statistical Power", fontweight="bold")
axes[1].set_ylabel("Mean Power")
axes[1].set_ylim(0, 0.60)

# Add numerical labels above each bar
for ax, bars in zip(axes, [bars1, bars2]):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{height:.3f}",
            (bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            fontsize=9
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle(
    "Consensus Voting Performance Across 200 Simulations",
    fontsize=13,
    fontweight="bold"
)

fig.tight_layout()

# Save high-resolution figures
fig.savefig(
    OUTPUT / "figure2_consensus_performance.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    OUTPUT / "figure2_consensus_performance.pdf",
    bbox_inches="tight"
)

plt.close(fig)

print("\nFigure 2 saved successfully!")
