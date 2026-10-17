from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = Path("erinmaley/data/consensus_report_7")
df = pd.read_csv(base / "consensus_signals_cdi_schubert.csv")

counts = df["SignificantMethodCount"].value_counts()
values = [int(counts.get(n, 0)) for n in [2, 3, 4]]
labels = ["2 of 4 methods", "3 of 4 methods", "4 of 4 methods"]
colors = ["#293B69", "#B28AA8", "#D49A76"]

assert sum(values) == len(df)

fig, ax = plt.subplots(figsize=(10, 3.2))
left = 0

for value, label, color in zip(values, labels, colors):
    ax.barh(0, value, left=left, color=color,
            height=0.5, label=label)
    ax.text(left + value / 2, 0, str(value),
            ha="center", va="center",
            color="white", fontsize=12, weight="bold")
    left += value

ax.set_xlim(0, len(df))
ax.set_yticks([])
ax.set_xlabel("Number of consensus-positive taxa")
ax.set_title(
    "Consensus Support Across 286 CDI Microbial Taxa",
    fontsize=14, weight="bold", pad=16
)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32),
          ncol=3, frameon=False)

for spine in ax.spines.values():
    spine.set_visible(False)

fig.tight_layout()
output = base / "figure3_cdi_consensus_distribution.png"
fig.savefig(output, dpi=300, bbox_inches="tight")
plt.close(fig)

print("CDI CONSENSUS FIGURE COMPLETE:", output)
print(dict(zip(labels, values)))
