"""
Functions for manipulating expression data.
"""

import sys
import numpy as np
import scipy.sparse
from skimage.filters import gaussian
from sklearn.linear_model import LinearRegression
import torch


def coo_to_torch(coo: scipy.sparse.coo_matrix, dense: bool = True) -> torch.Tensor:
    """Convert a scipy sparse COO matrix to a PyTorch tensor."""
    coo = np.ascontiguousarray(coo)
    indices = np.vstack((coo.row, coo.col))
    values = coo.data
    shape = coo.shape

    tensor = torch.sparse_coo_tensor(
        torch.LongTensor(indices), torch.FloatTensor(values), shape
    )

    return tensor.to_dense() if dense else tensor


def zero_pad(
    M: torch.Tensor, base_height: int, base_width: int, tile_size: int
) -> torch.Tensor:
    """Pad a matrix with zeros to ensure dimensions are multiples of `tile_size`."""
    horizontal_padding = (tile_size - base_width % tile_size) % tile_size
    vertical_padding = (tile_size - base_height % tile_size) % tile_size

    padded_width = base_width + horizontal_padding
    padded_height = base_height + vertical_padding

    padded = torch.zeros((padded_height, padded_width), dtype=M.dtype, device=M.device)
    padded[:base_height, :base_width] = M
    return padded


def tile_matrix(M: torch.Tensor, tile_size: int) -> torch.Tensor:
    """Break a matrix into tiles of size `tile_size`."""
    h = M.unfold(0, tile_size, tile_size)
    v = h.unfold(1, tile_size, tile_size)
    return v.permute(0, 2, 1, 3).reshape(-1, tile_size, tile_size)


def get_tile(matrix: torch.Tensor, index: int, tile_size: int) -> torch.Tensor:
    """Extract a tile from a matrix."""
    n = matrix.size(-1)
    row_start = (index * tile_size) // n * tile_size
    col_start = (index * tile_size) % n
    return matrix[
        ..., row_start : row_start + tile_size, col_start : col_start + tile_size
    ]


def get_landmarks_signal(
    level: int, tile_sizes: list, landmark_matrices: torch.Tensor, index: int = None
) -> torch.Tensor:
    """Calculate the landmark signal for a given level."""
    if level != 0 and index is None:
        raise ValueError("If level is not 0, index must be specified.")

    tile_size = tile_sizes[level]
    if level == 0:
        S = tile_matrix(landmark_matrices, tile_size).sum(dim=(-1, -2))
    else:
        M = get_tile(landmark_matrices, index, tile_sizes[level - 1])
        S = tile_matrix(M, tile_size).sum(dim=(-1, -2))

    return torch.nan_to_num(S / S.sum(dim=-1, keepdim=True))


def normalize_likelihood(
    likelihood: torch.Tensor, dim: int = 1, lambda_param: float = 0
) -> torch.Tensor:
    """Normalize the likelihood using sparsemax."""

    def normalize_column(col):
        return sparsemax(col, lambda_param=lambda_param)

    normalized = torch.stack(
        [normalize_column(likelihood[:, i]) for i in range(likelihood.size(dim))],
        dim=dim,
    )
    normalized[..., -1] = 0  # Set the last value to 0 for sparsity
    return normalized


def apply_mask(
    matrix: scipy.sparse.coo_matrix,
    mask: scipy.sparse.coo_matrix,
    return_dense: bool = False,
):
    """Apply a mask to a sparse matrix."""
    masked_matrix = scipy.sparse.coo_matrix(
        (matrix[mask.row, mask.col], (mask.row, mask.col)), shape=mask.shape
    )
    return masked_matrix.toarray() if return_dense else masked_matrix


def get_expectation(
    expression: torch.Tensor, probabilities: torch.Tensor
) -> torch.Tensor:
    """Calculate the expectation of a matrix based on probabilities."""
    return expression @ probabilities


def gaussian_filter(input_tensor: torch.Tensor, sigma: float = 1) -> torch.Tensor:
    """Apply a Gaussian filter to a tensor."""
    output = input_tensor.clone()
    for i in range(output.size(0)):
        max_val = input_tensor[i].max()
        filtered = gaussian(input_tensor[i].cpu().numpy(), sigma=sigma)
        output[i] = torch.tensor(filtered, device=input_tensor.device)
        if max_val > 0:
            output[i] *= max_val / output[i].max()
    return output


def get_threshold(input_tensor: torch.Tensor, k: int = 1024) -> float:
    """Calculate a threshold value based on an input tensor."""
    unique_vals = input_tensor[input_tensor > 0].unique()
    n = unique_vals.numel()
    m = n // 2

    model_1 = LinearRegression().fit(
        np.arange(m - k, m + k).reshape(-1, 1), unique_vals[m - k : m + k].numpy()
    )
    model_2 = LinearRegression().fit(
        np.arange(n - k, n).reshape(-1, 1), unique_vals[n - k :].numpy()
    )

    x_intersect = (model_2.intercept_ - model_1.intercept_) / (
        model_1.coef_ - model_2.coef_
    )
    return unique_vals[int(x_intersect)]


def sparsemax(z: torch.Tensor, lambda_param: float = 0) -> torch.Tensor:
    """Sparsemax function."""
    z_sorted, _ = torch.sort(z, descending=True)
    cumsum_z = torch.cumsum(z_sorted, dim=0)
    k_values = torch.arange(1, z.size(0) + 1, device=z.device)
    condition = 1 + k_values * z_sorted > cumsum_z
    k_z = k_values[condition][-1]
    tau_z = (cumsum_z[k_z - 1] - 1) / k_z
    return torch.clamp(z - tau_z, min=0)


def sparsegen(z: torch.Tensor, g=lambda x: x, lambda_param: float = 0) -> torch.Tensor:
    """Sparsegen function with custom generator."""
    g_z = g(z)
    lambda_min = 1 - torch.norm(g_z, p=1)
    lambda_param = max(lambda_param, lambda_min)
    return sparsemax(g_z / (1 - lambda_param))
