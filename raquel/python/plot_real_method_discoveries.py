
# BIOT 670I - Microbiome Abundance Consensus Machine
# Real-data Figure 1: Method Discoveries
# Dataset: cdi_vincent

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("results/real/cdi_vincent/consensus")
OUTPUT = Path("results/figures")
OUTPUT.mkdir(parents=True, exist_ok=True)

# Load both ANCOM-BC2 significance scenarios
robust = pd.read_csv(BASE / "method_discoveries_robust.csv")
adjusted = pd.read_csv(BASE / "method_discoveries_adjusted.csv")

methods = [
    "CLR_Wilcoxon",
    "ALDEx2",
    "ANCOMBC2",
    "DESeq2",
]

labels = [
    "CLR + Wilcoxon",
    "ALDEx2",
    "ANCOM-BC2",
    "DESeq2",
]

robust = robust.set_index("Method").loc[methods]
adjusted = adjusted.set_index("Method").loc[methods]

colors = ["#52796F", "#84A98C", "#D4A373", "#A98467"]

fig, ax = plt.subplots(figsize=(9, 5))

bars = ax.bar(
    labels,
    robust["Significant_OTUs"],
    color=colors,
    edgecolor="black",
    linewidth=0.5
)

# Show standard ANCOM-BC2 significance separately
ax.bar(
    "ANCOM-BC2",
    adjusted.loc["ANCOMBC2", "Significant_OTUs"],
    width=0.45,
    color="none",
    edgecolor="#8B4513",
    linewidth=1.5,
    linestyle="--",
    label="ANCOM-BC2 standard significance"
)

for bar in bars:
    height = bar.get_height()
    ax.annotate(
        f"{int(height)}",
        (bar.get_x() + bar.get_width() / 2, height),
        xytext=(0, 5),
        textcoords="offset points",
        ha="center",
        fontsize=10
    )

ax.annotate(
    "61 standard",
    (2, 61),
    xytext=(0, 5),
    textcoords="offset points",
    ha="center",
    fontsize=10
)

ax.set_title(
    "Differential Abundance Discoveries in the CDI Dataset",
    fontsize=13,
    fontweight="bold"
)

ax.set_ylabel("Number of Significant OTUs")
ax.set_ylim(0, 72)
ax.tick_params(axis="x", labelrotation=20)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(frameon=False)

fig.tight_layout()

fig.savefig(
    OUTPUT / "figure5_real_method_discoveries.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    OUTPUT / "figure5_real_method_discoveries.pdf",
    bbox_inches="tight"
)

plt.close(fig)

print("Real-data Figure 5 saved successfully!")

