
# BIOT 670I - Microbiome Abundance Consensus Machine
# Figure 7: Real-Data Method Agreement
# Dataset: cdi_vincent

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# Input and output directories
BASE = Path("results/real/cdi_vincent/consensus")
OUTPUT = Path("results/figures")
OUTPUT.mkdir(parents=True, exist_ok=True)

# Load primary consensus results
votes = pd.read_csv(BASE / "otu_votes_robust.csv")

# Keep OTUs supported by at least two methods
votes = votes.loc[votes["Votes"] >= 2].copy()

# Sort by agreement level, then OTU ID
votes = votes.sort_values(
    ["Votes", "Taxon"],
    ascending=[False, True]
)

# Methods used in the consensus
methods = [
    "CLR_Wilcoxon",
    "ALDEx2",
    "ANCOMBC2",
    "DESeq2"
]

method_labels = [
    "CLR + Wilcoxon",
    "ALDEx2",
    "ANCOM-BC2 (robust)",
    "DESeq2"
]

# Extract genus from taxonomy annotation
def get_genus(taxonomy):
    for part in str(taxonomy).split(";"):
        if part.startswith("g__"):
            genus = part.replace("g__", "")
            return genus if genus else "Unclassified"
    return "Unclassified"

votes["Genus"] = votes["Taxonomy"].apply(get_genus)

row_labels = [
    f"{row.Taxon} ({row.Genus})"
    for row in votes.itertuples()
]

# Convert method calls to numeric matrix
matrix = votes[methods].astype(int).to_numpy()

# Match existing capstone color palette
cmap = ListedColormap([
    "#E8E5DF",  # Not significant
    "#52796F"   # Significant
])

fig, ax = plt.subplots(figsize=(10, 5.5))

ax.imshow(
    matrix,
    cmap=cmap,
    vmin=0,
    vmax=1,
    aspect="auto",
    interpolation="nearest"
)

# Axis labels
ax.set_xticks(np.arange(len(methods)))
ax.set_xticklabels(method_labels, rotation=20, ha="right")

ax.set_yticks(np.arange(len(row_labels)))
ax.set_yticklabels(row_labels)

# Add checkmarks for significant findings
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        if matrix[i, j] == 1:
            ax.text(
                j, i, "✓",
                ha="center",
                va="center",
                fontsize=16,
                fontweight="bold",
                color="white"
            )

# Add subtle gridlines
ax.set_xticks(
    np.arange(-0.5, len(methods), 1),
    minor=True
)

ax.set_yticks(
    np.arange(-0.5, len(row_labels), 1),
    minor=True
)

ax.grid(
    which="minor",
    color="white",
    linewidth=2
)

ax.tick_params(which="minor", bottom=False, left=False)

# Legend
legend_items = [
    Patch(facecolor="#52796F", label="Significant"),
    Patch(facecolor="#E8E5DF", label="Not significant")
]

ax.legend(
    handles=legend_items,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.22),
    ncol=2,
    frameon=False
)

ax.set_title(
    "Agreement Across Consensus-Supported OTUs in the CDI Dataset",
    fontsize=13,
    fontweight="bold",
    pad=15
)

ax.set_xlabel("Differential Abundance Method")

fig.tight_layout()

# Save high-resolution figures
fig.savefig(
    OUTPUT / "figure7_real_method_agreement.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    OUTPUT / "figure7_real_method_agreement.pdf",
    bbox_inches="tight"
)

plt.close(fig)

print("Real-data Figure 7 saved successfully!")

