# BIOT 670I - Microbiome Abundance Consensus Machine
# Figure 1: DA Method Performance Across 200 Simulations

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Input and output directories
BASE = Path("results/da/combined_200")
OUTPUT = Path("results/figures")
OUTPUT.mkdir(parents=True, exist_ok=True)

# Differential abundance methods
METHODS = {
    "CLR + Wilcoxon": "clr_wilcoxon_performance.csv",
    "ALDEx2": "aldex2_performance.csv",
    "ANCOM-BC2": "ancombc_performance.csv",
    "DESeq2": "deseq2_performance.csv",
}

# Calculate mean performance across 200 simulations
records = []

for method, filename in METHODS.items():
    df = pd.read_csv(BASE / filename)

    records.append({
        "Method": method,
        "Mean_FDP": df["FDP"].mean(),
        "Mean_Power": df["Power"].mean(),
    })

summary = pd.DataFrame(records)

print("\nDA Method Performance (200 Simulations)")
print(summary.to_string(index=False))

# Create a two-panel comparison figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

colors = ["#52796F", "#84A98C", "#D4A373", "#A98467"]

# Panel A: Mean false discovery proportion
bars1 = axes[0].bar(
    summary["Method"],
    summary["Mean_FDP"],
    color=colors,
    edgecolor="black",
    linewidth=0.5
)

axes[0].set_title("A. False Discovery Proportion", fontweight="bold")
axes[0].set_ylabel("Mean FDP")
axes[0].set_ylim(0, 0.55)

# Panel B: Mean statistical power
bars2 = axes[1].bar(
    summary["Method"],
    summary["Mean_Power"],
    color=colors,
    edgecolor="black",
    linewidth=0.5
)

axes[1].set_title("B. Statistical Power", fontweight="bold")
axes[1].set_ylabel("Mean Power")
axes[1].set_ylim(0, 0.65)

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

    ax.tick_params(axis="x", labelrotation=25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle(
    "Differential Abundance Method Performance Across 200 Simulations",
    fontsize=13,
    fontweight="bold"
)

fig.tight_layout()

# Save high-resolution figures
fig.savefig(
    OUTPUT / "figure1_da_method_performance.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    OUTPUT / "figure1_da_method_performance.pdf",
    bbox_inches="tight"
)

plt.close(fig)

print("\nFigure 1 saved successfully!")
