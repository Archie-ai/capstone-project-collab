# BIOT 670I - Microbiome Abundance Consensus Machine
# Figure 3: False Discovery vs. Statistical Power

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Input and output directories
METHOD_BASE = Path("results/da/combined_200")
CONSENSUS_BASE = Path("results/da/consensus_v2_200")
OUTPUT = Path("results/figures")
OUTPUT.mkdir(parents=True, exist_ok=True)

# Individual differential abundance methods
METHODS = {
    "CLR + Wilcoxon": "clr_wilcoxon_performance.csv",
    "ALDEx2": "aldex2_performance.csv",
    "ANCOM-BC2": "ancombc_performance.csv",
    "DESeq2": "deseq2_performance.csv",
}

# Calculate mean FDP and power for individual methods
records = []

for method, filename in METHODS.items():
    df = pd.read_csv(METHOD_BASE / filename)

    records.append({
        "Approach": method,
        "Mean_FDP": df["FDP"].mean(),
        "Mean_Power": df["Power"].mean(),
        "Type": "Individual Method",
    })

# Calculate mean FDP and power for consensus voting rules
consensus = pd.read_csv(
    CONSENSUS_BASE / "consensus_performance_by_simulation.csv"
)

for rule in ["2-of-4", "3-of-4", "4-of-4"]:
    df = consensus[consensus["Rule"] == rule]

    records.append({
        "Approach": rule,
        "Mean_FDP": df["FDP"].mean(),
        "Mean_Power": df["Power"].mean(),
        "Type": "Consensus Rule",
    })

# Combine results
summary = pd.DataFrame(records)

print("\nFDP vs. Power Comparison (200 Simulations)")
print(summary.to_string(index=False))

# Create scatterplot
fig, ax = plt.subplots(figsize=(10, 6))

colors = {
    "Individual Method": "#52796F",
    "Consensus Rule": "#D4A373",
}

markers = {
    "Individual Method": "o",
    "Consensus Rule": "D",
}

for approach_type in ["Individual Method", "Consensus Rule"]:
    subset = summary[summary["Type"] == approach_type]

    ax.scatter(
        subset["Mean_FDP"],
        subset["Mean_Power"],
        color=colors[approach_type],
        marker=markers[approach_type],
        s=120,
        edgecolor="black",
        linewidth=0.7,
        label=approach_type,
        zorder=3,
    )

# Adjust label positions to prevent overlap
label_offsets = {
    "ALDEx2": (-10, 10),
    "2-of-4": (10, 10),
}

for _, row in summary.iterrows():
    offset = label_offsets.get(row["Approach"], (7, 7))

    ax.annotate(
        row["Approach"],
        (row["Mean_FDP"], row["Mean_Power"]),
        xytext=offset,
        textcoords="offset points",
        fontsize=9,
        ha="right" if row["Approach"] == "ALDEx2" else "left",
    )

# Format the figure
ax.set_title(
    "False Discovery vs. Statistical Power Across 200 Simulations",
    fontsize=13,
    fontweight="bold",
)

ax.set_xlabel("Mean False Discovery Proportion (FDP)")
ax.set_ylabel("Mean Statistical Power")

ax.set_xlim(-0.02, 0.50)
ax.set_ylim(0.25, 0.60)

ax.legend(frameon=False)
ax.grid(alpha=0.2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()

# Save high-resolution figures
fig.savefig(
    OUTPUT / "figure3_fdp_power_tradeoff.png",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    OUTPUT / "figure3_fdp_power_tradeoff.pdf",
    bbox_inches="tight"
)

plt.close(fig)

print("\nFigure 3 saved successfully!")
