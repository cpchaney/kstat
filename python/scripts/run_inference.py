"""
Run spatial inference by calculating bin-wise cell probability distributions.

Usage:
    python scripts/run_inference.py config/mouse_e18.5d.json
"""

import time

import pandas as pd
import torch
from kstat import model
from kstat.probability import compute_bin_inputs, load_config, prepare_inputs
from scipy import sparse
from scipy.io import mmwrite

# from kstat.utils import row_normalize, torch_mmread

# --------------------------------------------------
# Step 1: Configuration and device setup
# --------------------------------------------------

CONFIG_PATH = "../../config/mouse_e18.5d.json"  # Modify path as needed
K = 256  # Top-k expression values per gene
LAMBDA_PARAM = 0.1  # Regularization weight

# Select computation device
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device")

# --------------------------------------------------
# Step 2: Load configuration and data
# --------------------------------------------------

config = load_config(CONFIG_PATH)

# Load spatial landmark matrix, mask, and expression reference
landmark_counts, mask, ad, config = prepare_inputs(config, device)
landmarks = config["landmarks"]

# --------------------------------------------------
# Step 3: Compute bin-related inputs for modeling
# --------------------------------------------------

unique_bins, unique_bins_key, support, support_loc, spatial_expression = (
    compute_bin_inputs(
        landmark_counts=landmark_counts,
        mask=mask,
        landmarks=landmarks,
        common_landmarks=config["common_landmarks"],
        ad=ad,
        device=device,
    )
)

# --------------------------------------------------
# Step 4: Load reference expression matrix
# --------------------------------------------------

# Convert imputed expression matrix to torch tensor (genes x cells)
expression_tensor = torch.tensor(ad.layers["imputed"].T).float().to(device)

# For thresholding: compute per-gene K-th highest value
topk_values, _ = torch.topk(expression_tensor, K, dim=1)
thresholds = topk_values[:, -1]  # Currently unused, could filter low-expression genes

# --------------------------------------------------
# Step 5: Estimate parameters (mean and covariance)
# --------------------------------------------------

mean_values, covariances = model.infer_parameters_from_mixture(
    expression_tensor, k=None
)

# --------------------------------------------------
# Step 6: Compute bin assignment probabilities
# --------------------------------------------------

start_time = time.time()

cells_bins_probabilities = model.calculate_bin_probabilities(
    expression_data=expression_tensor,
    mean_values=mean_values,
    covariances=covariances,
    unique_bins=unique_bins,
    lambda_param=LAMBDA_PARAM,
)

execution_time = time.time() - start_time
print(f"Inference complete in {execution_time:.2f} seconds")

# --------------------------------------------------
# Step 7: Post-processing
# --------------------------------------------------

# Move result to CPU for saving or downstream use
cells_bins_probabilities = cells_bins_probabilities.cpu()

cells_bins_probabilities_sparse = sparse.coo_matrix(cells_bins_probabilities.numpy())

# Optionally: save to disk
# torch.save(cells_bins_probabilities, f"{config['project_root']}output/{config['experiment_name']}_cells_bins_probabilities.pt")
mmwrite(
    f"{config['project_root']}data/{config['experiment_name']}_cells_bins_probabilities.mtx",
    cells_bins_probabilities_sparse,
)

# Convert tensor to numpy, then to DataFrame
df = pd.DataFrame({"index": unique_bins_key.numpy()})

# Path to save the CSV file
output_path = (
    f"{config['project_root']}data/{config['experiment_name']}_unique_bins_key.csv"
)

# Write to CSV (without the default pandas index column)
df.to_csv(output_path, index=False)
