# BIOT 670I - Microbiome Abundance Consensus Machine
# Figure 4: Consensus Agreement Distribution Across 200 Simulations

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Input and output directories
BASE = Path("results/da/consensus_v2_200")
OUTPUT = Path("results/figures")
OUTPUT.mkdir(parents=True, exist_ok=True)

# Load consensus agreement results
agreement = pd.read_csv(
    BASE / "vote_distribution_by_simulation.csv"
)

# Calculate average number of taxa receiving each vote count
summary = (
    agreement.groupby("Votes")["Taxa_Count"]
    .mean()
    .reindex(range(5), fill_value=0)
    .reset_index()
)

print("\nConsensus Agreement Distribution (200 Simulations)")
print(summary.to_string(index=False))

# Create a two-panel consensus agreement figure
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

colors = [
    "#52796F",
    "#84A98C",
    "#B7C9A8",
    "#D4A373",
    "#A98467",
]

# Panel A: Complete distribution (0-4 votes)
bars1 = axes[0].bar(
    summary["Votes"].astype(str),
    summary["Taxa_Count"],
    color=colors,
    edgecolor="black",
    linewidth=0.5,
)

axes[0].set_title("A. Overall Agreement", fontweight="bold")
axes[0].set_ylabel("Mean Number of Taxa per Simulation")
axes[0].set_ylim(0, 105)

# Panel B: Taxa identified by at least one method
detected = summary[summary["Votes"] > 0]

bars2 = axes[1].bar(
    detected["Votes"].astype(str),
    detected["Taxa_Count"],
    color=colors[1:],
    edgecolor="black",
    linewidth=0.5,
)

axes[1].set_title("B. Agreement Among Detected Taxa", fontweight="bold")
axes[1].set_ylabel("Mean Number of Taxa per Simulation")
axes[1].set_ylim(0, 5)

# Add value labels and clean formatting
for ax, bars in zip(axes, [bars1, bars2]):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{height:.2f}",
            (bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            fontsize=9,
        )

    ax.set_xlabel("Number of Methods Identifying a Taxon")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle(
    "Differential Abundance Method Agreement Across 200 Simulations",
    fontsize=13,
    fontweight="bold",
)

fig.tight_layout()

# Save high-resolution figures
fig.savefig(
    OUTPUT / "figure4_consensus_agreement.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    OUTPUT / "figure4_consensus_agreement.pdf",
    bbox_inches="tight"
)

plt.close(fig)

print("\nFigure 4 saved successfully!")
