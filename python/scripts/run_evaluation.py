"""
Run evaluation of reconstructed gene expression using Sinkhorn divergence and cosine similarity.

Usage:
    python scripts/run_evaluation.py config/mouse_e18.5d.json
"""

import sys
import torch

from kstat.evaluation import (
    compute_cosine_similarity,
    compute_expected_expression,
    compute_sinkhorn_divergence,
)
from kstat.probability import (
    compute_bin_inputs,
    load_config,
    prepare_inputs,
)
from kstat.utils import row_normalize, torch_mmread

# --------------------------------------------------------
# Step 1: Load configuration and set device
# --------------------------------------------------------

# Accept config file path from command-line argument or use default
config_path = sys.argv[1] if len(sys.argv) > 1 else "config/mouse_e18.5d.json"
config = load_config(config_path)

# Use CUDA if available
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device")

# --------------------------------------------------------
# Step 2: Load input data (landmark matrix, mask, AnnData)
# --------------------------------------------------------

landmark_counts, mask, ad, config = prepare_inputs(config, device)
landmarks = config["landmarks"]

# Extract expression tensor (genes x cells)
expression_tensor = torch.tensor(ad.layers["imputed"].T).float().to(device)

# Mask landmark genes that are in the common set
masked_landmark_counts = landmark_counts.to_dense()[
    [gene in config["common_landmarks"] for gene in landmarks], :
] * mask.to_dense().reshape(-1)

# --------------------------------------------------------
# Step 3: Load model output: cells ? spatial bins probabilities
# --------------------------------------------------------

cells_bins_probabilities = torch.load(
    f"{config['project_root']}output/{config['experiment_name']}_cells_bins_probabilities.pt"
)

# Map each spatial bin to a unique binary pattern index
unique_bins_key = torch.unique(
    (masked_landmark_counts > 0).int(), dim=1, return_inverse=True
)[1]

# Create spatial support mask (1D) and support coordinates (2D)
support = mask.to_dense().reshape(-1) == 1
support = support.contiguous().to(device)
support_loc = mask.indices().float().T.contiguous().to(device)

# --------------------------------------------------------
# Step 4: Compute expected spatial expression from model
# --------------------------------------------------------

expected_expression = compute_expected_expression(
    expression_tensor,
    cells_bins_probabilities,
    unique_bins_key,
    support,
    device,
)

# Normalize the actual observed spatial expression
spatial_expression = row_normalize(masked_landmark_counts).to(device)

# --------------------------------------------------------
# Step 5: Evaluate reconstruction quality
# --------------------------------------------------------

sinkhorn_loss = compute_sinkhorn_divergence(
    expected_expression,
    spatial_expression[:, support],
    support_loc,
)

cosine_sim = compute_cosine_similarity(
    expected_expression,
    spatial_expression,
    support,
)

# --------------------------------------------------------
# Step 6: Output results
# --------------------------------------------------------

print("Mean Sinkhorn divergence:", sinkhorn_loss.mean().item())
print("Mean cosine similarity:", cosine_sim.mean().item())

# Save full vectors to disk
torch.save(
    sinkhorn_loss.cpu(),
    f"{config['project_root']}output/{config['experiment_name']}_sinkhorn.pt",
)
torch.save(
    cosine_sim.cpu(),
    f"{config['project_root']}output/{config['experiment_name']}_cosine.pt",
)
