import json

import numpy as np
import scanpy as sc
import torch

from kstat.utils import row_normalize, torch_mmread


def load_config(path):
    """
    Load a JSON configuration file.

    Args:
        path (str): Path to the JSON file.

    Returns:
        dict: Parsed configuration dictionary.
    """
    with open(path, "r") as file:
        return json.load(file)


def prepare_inputs(config, device):
    """
    Load and prepare key input data for modeling.

    This includes:
      - Reading the spatial landmark matrix and binary mask.
      - Loading the single-nucleus RNA-seq reference data.
      - Filtering for common genes between spatial and reference datasets.

    Args:
        config (dict): Configuration dictionary.
        device (torch.device): Target device for tensors.

    Returns:
        tuple: (landmark_counts, mask, adata, updated_config)
    """
    # Load spatial landmark count matrix
    landmark_counts = torch_mmread(
        f"{config['project_root']}data/{config['experiment_name']}_landmarks_matrix.mtx"
    ).to_sparse_csr()

    # Load binary mask (scaled to spatial resolution)
    mask = torch_mmread(
        f"{config['project_root']}data/{config['experiment_name']}_scaled_mask.mtx"
    ).coalesce()

    # Load reference expression matrix (imputed)
    ad = sc.read_h5ad(
        f"{config['project_root']}../kidney_development/data/{config['experiment_name']}.h5ad"
    )

    # Keep only genes with nonzero expression across all cells
    nonzero_genes = np.array(ad.layers["imputed"].sum(axis=0) > 0).flatten()
    ad = ad[:, nonzero_genes]

    # Keep only landmark genes present in the reference dataset
    landmarks = config["landmarks"]
    common_landmarks = [gene for gene in landmarks if gene in ad.var_names]
    config["common_landmarks"] = common_landmarks
    ad = ad[:, common_landmarks]

    return landmark_counts, mask, ad, config


def compute_bin_inputs(landmark_counts, mask, landmarks, common_landmarks, ad, device):
    """
    Process spatial measurements into binary landmark presence/absence patterns (bins),
    and extract associated spatial information (mask, support coordinates).

    Args:
        landmark_counts (torch.Tensor): Landmark count matrix [genes x pixels].
        mask (torch.Tensor): Binary mask indicating valid spatial bins.
        landmarks (list): Full list of configured landmark genes.
        common_landmarks (list): Landmark genes present in both spatial and reference data.
        ad (AnnData): Filtered reference expression matrix.
        device (torch.device): Target device (e.g., CUDA or CPU).

    Returns:
        tuple:
            - unique_bins (Tensor): Unique binary presence/absence patterns.
            - support (Tensor): Binary tensor indicating active spatial locations.
            - support_loc (Tensor): 2D spatial coordinates for each active bin.
            - spatial_expression (Tensor): Normalized spatial signal for each landmark.
    """
    # Restrict spatial landmark counts to shared genes
    landmark_mask = [gene in common_landmarks for gene in landmarks]
    masked_landmark_counts = landmark_counts.to_dense()[
        landmark_mask, :
    ] * mask.to_dense().reshape(-1)

    # Convert counts to binary indicators (presence/absence)
    binarized_landmarks = (masked_landmark_counts > 0).to(torch.int)

    # Identify all unique presence/absence patterns across spatial bins
    unique_bins, _ = torch.unique(binarized_landmarks, dim=1, return_inverse=True)

    # Create binary mask of active bins and convert to target device
    support = (mask.to_dense().reshape(-1) == 1).contiguous().to(device)

    # Get spatial coordinates of nonzero bins
    support_loc = mask.indices().float().T.contiguous().to(device)

    # Extract landmark signals at supported locations and normalize
    spatial_expression = row_normalize(binarized_landmarks[:, support.cpu()]).float()
    spatial_expression = spatial_expression.contiguous().to(device)

    return unique_bins.to(device), support, support_loc, spatial_expression
