"""
Spatially-Aware Gene Expression Modeling Module

This module provides functionality to model and infer spatial patterns in gene
expression data using a combination of Gaussian Mixture Models (GMM), probabilistic
inference, and graph-based smoothing techniques.

It includes methods for:
- Fitting GMMs to expression profiles and selecting optimal components using BIC.
- Estimating statistical parameters (means and covariances) for landmark gene expression.
- Calculating probabilistic assignments of bins (spatial units) to inferred states.
- Applying spatial and cell-wise smoothing to probability distributions.
- Defining a Bayesian spatial model using Pyro for probabilistic inference.
- Providing a guide function for variational inference with Pyro.

Functions:
- `fit_gmm_with_best_bic`: Fits GMMs with varying components and selects the best one using BIC.
- `compute_covariance_matrix`: Computes the empirical covariance matrix of expression data.
- `infer_parameters_from_mixture`: Infers landmark-specific mean and covariance values using GMMs.
- `calculate_bin_probability`: Calculates the likelihood of a spatial bin being in a specific expression state.
- `calculate_bin_probabilities`: Computes probabilities for all bins given expression and landmark statistics.
- `smooth_probability_torch`: Smooths probability estimates using adjacency graphs (spatial and cell-cell).
- `bayesian_spatial_model`: Defines a generative probabilistic model for spatially-aware expression data.
- `guide`: Variational guide for Bayesian inference, modeling latent expression probabilities.

Dependencies:
- numpy
- torch
- pyro
- scikit-learn (GaussianMixture)
- tqdm
- statlas.utils (must define `sparsegen` for soft assignment)

This module is designed to support downstream spatial transcriptomics analyses,
where accurate probabilistic modeling of gene expression across spatial contexts
can reveal biologically meaningful patterns.

Example use:
    1. Estimate mean and covariance from expression data.
    2. Use bin probabilities to assign spatial states.
    3. Apply spatial smoothing.
    4. Perform Bayesian inference with Pyro to model posterior distributions.

Author: Christopher Chaney
"""

import numpy as np
import pyro
import pyro.distributions as dist
import torch
import torch.nn.functional as F
import tqdm
from pyro.distributions import Normal, constraints
from sklearn.mixture import GaussianMixture

import statlas.utils


def fit_gmm_with_best_bic(data, min_components=2, max_components=10):
    best_bic = float("inf")
    best_model = None
    bic_scores = []

    for k in range(min_components, max_components + 1):
        # Fit a Gaussian Mixture Model with k components
        gmm = GaussianMixture(n_components=k, random_state=0)
        gmm.fit(data)

        # Calculate BIC
        bic = gmm.bic(data)
        bic_scores.append(bic)
        print(f"Components: {k}, BIC: {bic}")

        # Update the best model if the BIC is lower
        if bic < best_bic:
            best_bic = bic
            best_model = gmm

    return best_model, best_bic, bic_scores


def compute_covariance_matrix(X):
    """
    Computes the covariance matrix of a dataset using PyTorch.

    Args:
    X (torch.Tensor): A tensor of shape (n, d), where n is the number of samples and d is the number of features.

    Returns:
    torch.Tensor: Covariance matrix of shape (d, d).
    """
    n_samples = X.shape[0]

    # Compute the mean of each column (feature-wise mean)
    mean_X = torch.mean(X, dim=0, keepdim=True)

    # Center the data
    X_centered = X - mean_X

    # Compute covariance matrix
    cov_matrix = (X_centered.T @ X_centered) / (n_samples - 1)

    return cov_matrix


def infer_parameters_from_mixture(expression_tensor, k=2):
    """Infer parameters for the given expression data and landmarks.

    Args:
        expression_tensor (torch.Tensor): Expression data indexed by landmark
            names, with columns representing gene expressions.
    """

    num_landmarks = expression_tensor.shape[0]
    device = expression_tensor.device

    # Instead of in-place modifications, we will build the tensors dynamically.
    means = []
    covariances = []

    for i in range(num_landmarks):
        x = expression_tensor[i, :].cpu().numpy()

        if k is not None:
            gmm = GaussianMixture(n_components=k, random_state=42)
            gmm.fit(x.reshape(-1, 1))
        else:
            gmm, best_bic, bic_scores = fit_gmm_with_best_bic(
                x.reshape(-1, 1), max_components=8
            )

        components = gmm.predict(x.reshape(-1, 1))
        on_component = np.argmax(gmm.means_)

        means.append(
            torch.tensor(
                np.insert(
                    gmm.means_[on_component], 0, x[components != on_component].mean()
                )
            )
        )
        covariances.append(
            torch.tensor(
                np.insert(
                    gmm.covariances_[on_component],
                    0,
                    x[components != on_component].std() ** 2,
                )
            )
        )

    # Stack the lists into tensors
    means = torch.stack(means).squeeze().to(device)
    covariances = torch.stack(covariances).squeeze().to(device)

    return means, covariances


def calculate_bin_probability(
    expression_data,
    mean_values,
    covariances,
    # correlation_matrix,
    landmark_indices,
    lambda_param,
):
    """Get the probability of a bin being in a particular state.

    Args:
        expression_data (torch.Tensor): Expression data.
        mean_values (torch.Tensor): Mean values for each landmark.
        covariancess (torch.Tensor): Covariances for each landmark.
        landmark_indices (torch.Tensor): Landmark indices.

    Returns:
        torch.Tensor: Probability of the bin being in a particular state.
    """

    mean_selected = mean_values.gather(1, landmark_indices.long()[:, None]).reshape(-1)
    combined_covariance = torch.diag(
        covariances.gather(1, landmark_indices.long()[:, None]).reshape(-1)
    )
    # covariances_selected = covariances.gather(
    #     1, landmark_indices.long()[:, None]
    # ).reshape(-1)
    # combined_covariance = (
    #     torch.outer(covariances_selected, covariances_selected) * correlation_matrix
    # )

    mvn = torch.distributions.MultivariateNormal(
        loc=mean_selected, covariance_matrix=combined_covariance
    )
    log_probs = mvn.log_prob(expression_data.T)

    maximum_magnitude = np.floor(np.log(torch.finfo(log_probs.dtype).max))
    clipped_log_probs = torch.clamp(
        log_probs, min=-maximum_magnitude, max=maximum_magnitude
    )

    # jacobian = torch.prod(1 / (expression_data + 1), dim=0)
    probabilities = torch.exp(clipped_log_probs)
    # probabilities = probabilities * jacobian

    if probabilities.sum() == 0:
        probabilities = torch.zeros(log_probs.shape)
    else:
        probabilities = probabilities / probabilities.sum()
        probabilities = statlas.utils.sparsegen(
            probabilities, lambda_param=lambda_param
        )
        # probabilities = probabilities / probabilities.sum()

    return probabilities


def calculate_bin_probabilities(
    expression_data,
    mean_values,
    covariances,
    unique_bins,
    lambda_param,
):
    """Calculate bin probabilities for each unique bin.

    Args:
        expression_data (torch.Tensor): Expression data.
        mean_values (torch.Tensor): Mean values for each landmark.
        covariances (torch.Tensor): Covariances matrices for each landmark.
        unique_bins (torch.Tensor): Unique bin indices.

    Returns:
        torch.Tensor: Bin probabilities.
    """
    device = expression_data.device
    mean_values.to(device=device)
    covariances.to(device=device)
    unique_bins.to(device=device)

    # covariance_matrix = compute_covariance_matrix(expression_data.T)
    # std_dev = torch.sqrt(torch.diag(covariance_matrix))

    # correlation_matrix = covariance_matrix / torch.outer(std_dev, std_dev)

    bin_probabilities = torch.zeros(expression_data.shape[1], unique_bins.shape[1]).to(
        device
    )

    for j in tqdm.tqdm(range(unique_bins.shape[1]), position=0, leave=True):
        landmark_indices = unique_bins[:, j].long()
        bin_probabilities[:, j] = calculate_bin_probability(
            expression_data,
            mean_values,
            covariances,
            # correlation_matrix,
            landmark_indices,
            lambda_param,
        )

    bin_probabilities[:, 0] = 0

    return bin_probabilities


def smooth_probability_torch(
    P_initial, spatial_adj, cell_adj, lambda_spatial=0.5, lambda_cell=0.5, num_iters=20
):
    """
    Smooths probability estimates using spatial and cell-cell graph adjacency matrices.

    Parameters:
    - P_initial: (N_cells, N_bins) torch.Tensor, initial probability estimates.
    - spatial_adj: (N_bins, N_bins) torch.sparse.Tensor, adjacency matrix for bin neighbors.
    - cell_adj: (N_cells, N_cells) torch.sparse.Tensor, adjacency matrix for nearest-neighbor cells.
    - lambda_spatial: float, weight for spatial smoothing.
    - lambda_cell: float, weight for cell-cell consistency.
    - num_iters: int, number of iterations for smoothing.

    Returns:
    - P_smooth: (N_cells, N_bins) torch.Tensor, smoothed probability estimates.
    """
    P_smooth = P_initial.clone()

    for _ in range(num_iters):
        # Spatial smoothing: propagate over neighboring bins
        P_spatial = torch.sparse.mm(spatial_adj, P_smooth)
        P_spatial /= torch.clamp(spatial_adj.sum(dim=1, keepdim=True), min=1e-6)

        # Cell-cell smoothing: propagate over nearest-neighbor cells
        P_neighbors = torch.sparse.mm(cell_adj, P_smooth)
        P_neighbors /= torch.clamp(cell_adj.sum(dim=1, keepdim=True), min=1e-6)

        # Weighted combination
        P_smooth = (
            (1 - lambda_spatial - lambda_cell) * P_initial
            + lambda_spatial * P_spatial
            + lambda_cell * P_neighbors
        )

    return P_smooth


def bayesian_spatial_model_batched(
    P_initial_batch,
    E_sn_batch,
    cell_adj_batch,
    norm_cell_batch,
    spatial_adj,
    norm_spatial,
    batch_idx,
    *,
    lambda_data=1.0,
    lambda_spatial=1.0,
    lambda_cell=1.0,
):
    device = P_initial_batch.device
    batch_size, N_bins = P_initial_batch.shape

    device = P_initial_batch.device

    with pyro.plate("cells", size=batch_size):
        logits = pyro.sample(
            "logits",
            dist.Normal(
                torch.tensor(0.0, device=device), torch.tensor(1.0, device=device)
            )
            .expand([N_bins])
            .to_event(1),
        )

    # Convert logits to probability distribution over bins
    P_latent = F.softmax(logits, dim=-1)
    P_latent = torch.clamp(P_latent, min=1e-6)  # avoid numerical instability

    # --- Reconstruct spatial gene expression ---
    E_pred = torch.matmul(P_latent.T, E_sn_batch)  # [N_bins x N_genes]
    E_smooth = torch.matmul(spatial_adj, E_pred)
    spatial_loss = torch.sum((E_pred - norm_spatial.unsqueeze(1) * E_smooth) ** 2)

    # --- Cell-cell smoothing ---
    P_neighbors = torch.sparse.mm(cell_adj_batch, P_latent)
    cell_loss = torch.sum((P_latent - norm_cell_batch.unsqueeze(1) * P_neighbors) ** 2)

    # --- Fit to initial probability estimate ---
    data_fit_loss = torch.sum((P_latent - P_initial_batch) ** 2)

    total_loss = (
        lambda_data * data_fit_loss
        + lambda_spatial * spatial_loss
        + lambda_cell * cell_loss
    )

    pyro.factor("loss", -total_loss)  # ELBO-friendly loss site


def guide_batched(
    P_initial_batch,
    E_sn_batch,
    cell_adj_batch,
    norm_cell_batch,
    spatial_adj,
    norm_spatial,
    batch_idx,
    *,
    lambda_data=1.0,
    lambda_spatial=1.0,
    lambda_cell=1.0,
):
    device = P_initial_batch.device
    batch_size, N_bins = P_initial_batch.shape

    # Variational parameters (unconstrained)
    loc_full = pyro.param("logits_loc")  # shape [N_cells, N_bins]
    scale_full = pyro.param(
        "logits_scale",
        torch.ones_like(loc_full) * 0.1,
        constraint=constraints.positive,
    )

    loc = loc_full[batch_idx].to(device)
    scale = scale_full[batch_idx].to(device)

    with pyro.plate("cells", size=batch_size):
        logits = pyro.sample(
            "logits", Normal(loc, scale).to_event(1)
        )  # shape [batch_size, N_bins]
