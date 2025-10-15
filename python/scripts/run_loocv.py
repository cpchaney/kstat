"""
Run leave-one-out cross-validation (LOOCV) to assess spatial reconstruction accuracy.

Usage:
    python scripts/run_loocv.py config/mouse_e18.5d.json
"""

import sys
import torch

from kstat.evaluation import run_loocv_evaluation, save_loocv_outputs
from kstat.probability import (
    compute_bin_inputs,
    load_config,
    prepare_inputs,
)
from kstat.utils import row_normalize, torch_mmread

# -----------------------------------------------------
# Step 1: Load configuration and set up device
# -----------------------------------------------------

# Allow optional config path via CLI
config_path = sys.argv[1] if len(sys.argv) > 1 else "config/mouse_e18.5d.json"
config = load_config(config_path)

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device")

# -----------------------------------------------------
# Step 2: Load spatial and reference expression data
# -----------------------------------------------------

# landmark_counts: landmark matrix (genes × pixels)
# mask: binary spatial mask
# ad: AnnData reference object
landmark_counts, mask, ad, config = prepare_inputs(config, device)

landmarks = config["landmarks"]
expression_tensor = torch.tensor(ad.layers["imputed"].T).float().to(device)

# -----------------------------------------------------
# Step 3: Compute bin inputs
# -----------------------------------------------------

unique_bins, support, support_loc, spatial_expression = compute_bin_inputs(
    landmark_counts,
    mask,
    landmarks,
    config["common_landmarks"],
    ad,
    device,
)

# -----------------------------------------------------
# Step 4: Load model output (cell ? bin probabilities)
# -----------------------------------------------------

cells_bins_probabilities = torch.load(
    f"{config['project_root']}output/{config['experiment_name']}_cells_bins_probabilities.pt"
)

# Map each pixel to a unique bin index
unique_bins_key = torch.unique(
    (landmark_counts.to_dense() > 0).int(), dim=1, return_inverse=True
)[1]

# -----------------------------------------------------
# Step 5: Run LOOCV evaluation
# -----------------------------------------------------

results_df = run_loocv_evaluation(
    ad=ad,
    common_landmarks=config["common_landmarks"],
    landmark_counts=landmark_counts,
    mask=mask,
    expression_tensor=expression_tensor,
    config=config,
    LAMBDA_PARAM=config.get("lambda", 0.1),
    support=support,
    support_loc=support_loc,
    spatial_expression=spatial_expression,
)

# Save per-gene evaluation results
results_df.to_csv(
    f"{config['project_root']}data/{config['experiment_name']}_loocv_results.csv"
)

# -----------------------------------------------------
# Step 6: Save LOOCV matrices for worst-performing genes
# -----------------------------------------------------

# Get top 10 highest Sinkhorn divergence genes (worst reconstructions)
error_genes = (
    results_df.sort_values("sinkhorn_divergence", ascending=False)
    .head(10)
    .index.tolist()
)

save_loocv_outputs(
    ad=ad,
    config=config,
    error_landmarks=error_genes,
    landmark_counts=landmark_counts,
    mask=mask,
    expression_tensor=expression_tensor,
    support_loc=support_loc,
    spatial_expression=spatial_expression,
    cells_bins_probabilities=cells_bins_probabilities,
    unique_bins_key=unique_bins_key,
)
