
# BIOT 670I - Microbiome Abundance Consensus Machine
# Figure 6: Real-Data Consensus Results
# Dataset: cdi_vincent

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Input and output directories
BASE = Path("results/real/cdi_vincent/consensus")
OUTPUT = Path("results/figures")
OUTPUT.mkdir(parents=True, exist_ok=True)

# Load consensus results
robust = pd.read_csv(
    BASE / "consensus_summary_robust.csv"
)

adjusted = pd.read_csv(
    BASE / "consensus_summary_adjusted.csv"
)

# Ensure consistent threshold order
rules = ["2-of-4", "3-of-4", "4-of-4"]

robust = robust.set_index("Rule").loc[rules]
adjusted = adjusted.set_index("Rule").loc[rules]

# Prepare bar positions
x = np.arange(len(rules))
width = 0.35

# Match existing capstone figure colors
colors = {
    "robust": "#52796F",
    "adjusted": "#D4A373"
}

fig, ax = plt.subplots(figsize=(9, 5))

bars1 = ax.bar(
    x - width / 2,
    robust["Significant_OTUs"],
    width,
    label="ANCOM-BC2 robust",
    color=colors["robust"],
    edgecolor="black",
    linewidth=0.5
)

bars2 = ax.bar(
    x + width / 2,
    adjusted["Significant_OTUs"],
    width,
    label="ANCOM-BC2 standard",
    color=colors["adjusted"],
    edgecolor="black",
    linewidth=0.5
)

# Add numerical labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()

        ax.annotate(
            f"{int(height)}",
            (
                bar.get_x() + bar.get_width() / 2,
                height
            ),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            fontsize=10
        )

ax.set_title(
    "Consensus Discoveries in the CDI Dataset",
    fontsize=13,
    fontweight="bold"
)

ax.set_xlabel("Consensus Threshold")
ax.set_ylabel("Number of Consensus-Supported OTUs")

ax.set_xticks(x)
ax.set_xticklabels(rules)
ax.set_ylim(0, 30)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(frameon=False)

fig.tight_layout()

# Save publication-quality figures
fig.savefig(
    OUTPUT / "figure6_real_consensus.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    OUTPUT / "figure6_real_consensus.pdf",
    bbox_inches="tight"
)

plt.close(fig)

print("Real-data Figure 6 saved successfully!")

