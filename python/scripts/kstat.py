"""
KSTAT algorithm for kidney CARTANA ISS and sn-RNA-seq data processing.
"""

import os
import json
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import torch
import tqdm
import matplotlib.pyplot as plt
import geomloss
from scipy.io import mmwrite
from typing import Tuple

import kstat.model
import kstat.utils


def torch_mmread(file_path: str) -> torch.Tensor:
    """Read a Matrix Market file and convert it to a PyTorch sparse tensor."""
    coo = scipy.io.mmread(file_path)
    indices = np.vstack((coo.row, coo.col))
    values = coo.data
    return torch.sparse_coo_tensor(
        torch.tensor(indices, dtype=torch.int64),
        torch.tensor(values, dtype=torch.float32),
        coo.shape,
    )


def row_normalize(tensor: torch.Tensor) -> torch.Tensor:
    """Normalize rows of a tensor."""
    row_sums = tensor.sum(dim=1, keepdim=True)
    row_sums[row_sums == 0] = 1.0
    return tensor / row_sums


def load_config(config_path: str) -> dict:
    """Load configuration file."""
    with open(config_path, "r") as file:
        return json.load(file)


def save_config(config: dict, config_path: str):
    """Save configuration file."""
    with open(config_path, "w") as file:
        json.dump(config, file, indent=4)


def filter_expression_data(ad, landmarks: list) -> Tuple[torch.Tensor, list]:
    """Filter expression data for common landmarks."""
    common_landmarks = [landmark for landmark in landmarks if landmark in ad.var_names]
    ad = ad[:, common_landmarks]
    expression_tensor = torch.tensor(ad.layers["imputed"]).float()
    return expression_tensor, common_landmarks


def compute_bin_probabilities(
    expression_tensor, mean_values, covariances, unique_bins, lambda_param
):
    """Compute bin probabilities for cells."""
    return kstat.model.calculate_bin_probabilities(
        expression_tensor,
        mean_values,
        covariances,
        unique_bins,
        lambda_param,
    )


def compute_sinkhorn_loss(
    expected_expression: torch.Tensor,
    spatial_expression: torch.Tensor,
    support_loc: torch.Tensor,
    device: str,
) -> torch.Tensor:
    """Calculate Sinkhorn divergence loss."""
    loss_fn = geomloss.SamplesLoss(loss="sinkhorn", p=2, blur=0.05, scaling=0.7)
    loss = torch.zeros(expected_expression.shape[0], device=device)

    for i in tqdm.tqdm(range(expected_expression.shape[0])):
        loss[i] = loss_fn(
            expected_expression[i, :],
            support_loc,
            spatial_expression[i, :],
            support_loc,
        )
    return loss


def compute_similarity(
    expected_expression: torch.Tensor,
    spatial_expression: torch.Tensor,
    support: torch.Tensor,
) -> torch.Tensor:
    """Calculate cosine similarity."""
    similarity = torch.zeros(expected_expression.shape[0])
    for i in tqdm.tqdm(range(expected_expression.shape[0])):
        similarity[i] = torch.nn.functional.cosine_similarity(
            expected_expression[i, :][None, :],
            spatial_expression[i, support][None, :],
        )
    return similarity


def save_results(
    config: dict,
    cells_bins_probabilities: torch.Tensor,
    unique_bins_key: torch.Tensor,
    unique_bins: torch.Tensor,
    loss: torch.Tensor,
    common_landmarks: list,
):
    """Save processed data to files."""
    # Save cells-bins probabilities
    mmwrite(
        f"{config['project_root']}data/{config['experiment_name']}_cells_bins_probabilities.mtx",
        scipy.sparse.coo_matrix(cells_bins_probabilities.cpu().numpy()),
    )

    # Save unique bins and keys
    pd.DataFrame(unique_bins_key.cpu().numpy(), columns=["index"]).to_csv(
        f"{config['project_root']}data/{config['experiment_name']}_unique_bins_key.csv",
        index=False,
    )
    mmwrite(
        f"{config['project_root']}data/{config['experiment_name']}_unique_bins.mtx",
        scipy.sparse.coo_matrix(unique_bins.cpu().numpy()),
    )

    # Save Wasserstein loss
    pd.DataFrame({"gene": common_landmarks, "loss": loss.cpu().numpy()}).to_csv(
        f"{config['project_root']}data/mouse_e15.5d_kidney_wasserstein_loss.csv",
        index=False,
    )


def main(config_path: str):
    config = load_config(config_path)

    # Set device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using {device} device")
    torch.set_default_dtype(torch.float32)

    # Load data
    landmark_counts = torch_mmread(
        f"{config['project_root']}data/{config['experiment_name']}_landmarks_matrix.mtx"
    ).to_sparse_csr()
    mask = torch_mmread(
        f"{config['project_root']}data/{config['experiment_name']}_scaled_mask.mtx"
    ).coalesce()
    ad = sc.read_h5ad(
        f"{config['project_root']}../kidney_development/data/{config['experiment_name']}.h5ad"
    )

    # Process data
    expression_tensor, common_landmarks = filter_expression_data(
        ad, config["landmarks"]
    )
    config["common_landmarks"] = common_landmarks

    masked_landmark_counts = landmark_counts.to_dense()[
        [landmark in common_landmarks for landmark in config["landmarks"]], :
    ] * mask.to_dense().reshape(-1)
    binarized_landmark_counts = (masked_landmark_counts > 0).to(torch.int)

    unique_bins, unique_bins_key = torch.unique(
        binarized_landmark_counts, dim=1, return_inverse=True
    )
    support = mask.to_dense().reshape(-1).bool().to(device)
    support_loc = mask.indices().float().T.contiguous().to(device)

    # Infer parameters
    mean_values, covariances = kstat.model.infer_parameters_from_mixture(
        expression_tensor, k=None
    )

    # Compute bin probabilities
    start_time = time.time()
    cells_bins_probabilities = compute_bin_probabilities(
        expression_tensor, mean_values, covariances, unique_bins.to(device), 0
    )
    print(f"Execution time: {time.time() - start_time} seconds")

    # Calculate Sinkhorn loss
    expected_expression = expression_tensor @ cells_bins_probabilities.cpu()
    expected_expression = row_normalize(expected_expression[:, unique_bins_key.cpu()])
    loss = compute_sinkhorn_loss(
        expected_expression.to(device),
        row_normalize(masked_landmark_counts).to(device),
        support_loc,
        device,
    )

    # Save results
    save_results(
        config,
        cells_bins_probabilities,
        unique_bins_key,
        unique_bins,
        loss,
        common_landmarks,
    )
    save_config(config, config_path)


if __name__ == "__main__":
    CONFIG_PATH = "/path/to/config.json"
    main(CONFIG_PATH)
