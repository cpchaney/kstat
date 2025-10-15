# algorithm_illustration.py

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity

# ------------------------------------------------------------
# STEP 1: Define spatial bins and simulated cell expression
# ------------------------------------------------------------

# Define binary landmark expression signatures for each spatial bin
bins = pd.DataFrame(
    {
        "Bin1": [1, 0, 0],  # expresses GeneA only
        "Bin2": [0, 1, 0],  # expresses GeneB only
        "Bin3": [0, 0, 1],  # expresses GeneC only
        "Bin4": [0, 1, 1],  # expresses GeneB and GeneC
    },
    index=["GeneA", "GeneB", "GeneC"]
)

# Define normalized expression levels for 2 cells across the same genes
cells = pd.DataFrame(
    {
        "CellX": [0.9, 0.2, 0.1],  # primarily GeneA
        "CellY": [0.1, 0.6, 0.8],  # primarily GeneB & GeneC
    },
    index=["GeneA", "GeneB", "GeneC"]
)

# ------------------------------------------------------------
# STEP 2: Compute cosine similarity and bin assignment probabilities
# ------------------------------------------------------------

# Cosine similarity: cell profiles vs. bin landmark patterns
bin_profiles = bins.T.values   # shape: (4 bins, 3 genes)
cell_profiles = cells.T.values  # shape: (2 cells, 3 genes)

# Compute cosine similarity between each cell and each bin
similarity = cosine_similarity(cell_profiles, bin_profiles)

# Normalize row-wise to get probability distribution per cell across bins
probabilities = similarity / similarity.sum(axis=1, keepdims=True)

# Wrap in DataFrame for inspection
prob_df = pd.DataFrame(probabilities, index=cells.columns, columns=bins.columns)

# Optional: normalize again column-wise for visualization (across cells)
prob_df_norm = prob_df.div(prob_df.sum(axis=0), axis=1)

# ------------------------------------------------------------
# STEP 3: Heatmap visualization
# ------------------------------------------------------------

plt.figure(figsize=(6, 3))
sns.heatmap(prob_df_norm, annot=True, cmap="YlGnBu", cbar=True)
plt.title("Normalized Probability of Cell Origin by Bin")
plt.xlabel("Spatial Bins")
plt.ylabel("Cells")
plt.tight_layout()
plt.savefig("../../image/cell_bin_heatmap.png", dpi=300)
plt.show()

# ------------------------------------------------------------
# STEP 4: Bar plot visualization of probabilities
# ------------------------------------------------------------

# Transpose for plotting: bins on x-axis
ax = prob_df_norm.T.plot(kind="bar", figsize=(6, 4), alpha=0.9)
plt.title("Normalized Bin Probabilities per Cell")
plt.ylabel("Probability")
plt.xlabel("Spatial Bins")
plt.ylim(0, 1.0)
plt.legend(title="Cells")
plt.tight_layout()
plt.savefig("../../image/cell_bin_barplot.png", dpi=300)
plt.show()

# ------------------------------------------------------------
# STEP 5: Simulated expression bar plots per cell
# ------------------------------------------------------------

# Assign colors to genes
gene_colors = {
    "Gene A": "#C10000",  # Red
    "Gene B": "#09B250",  # Green
    "Gene C": "#0170C0",  # Blue
    "Gene D": "#FFFE00",  # Yellow
}

# Simulated gene expression per cell
cell_expressions = {
    "Cell X": {"Gene A": 0.9, "Gene B": 0.2, "Gene C": 0.1, "Gene D": 0.6},
    "Cell Y": {"Gene A": 0.1, "Gene B": 0.6, "Gene C": 0.8, "Gene D": 0.3},
}

# Bar chart visual settings
bar_width = 0.8
spacing = 1.2
height = 3

# Generate bar plots for each cell
for cell_name, gene_data in cell_expressions.items():
    genes = list(gene_data.keys())
    values = list(gene_data.values())
    colors = [gene_colors[gene] for gene in genes]

    n_genes = len(genes)
    width = n_genes * spacing

    plt.figure(figsize=(width, height))
    bars = plt.bar(genes, values, color=colors, width=bar_width)

    # Set aesthetics
    plt.ylim(0, 1)
    plt.ylabel("Expression")
    plt.xticks(rotation=0)
    plt.tight_layout()

    # Save with transparent background
    filename = f"../../image/{cell_name.replace(' ', '_').lower()}_expression.png"
    plt.savefig(filename, dpi=300, transparent=True)
    plt.close()

    print(f"Saved: {filename}")
