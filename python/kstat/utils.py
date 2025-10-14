"""
Utilities for manipulating and analyzing transcriptomics data using PyTorch and SciPy.
"""

import sys

import numpy as np
import scipy.sparse
import skimage.filters
import sklearn.linear_model
import torch
import torch.nn.functional as F
from scipy.io import mmread, mmwrite


def torch_mmread(file_path):
    """
    Load a Matrix Market (.mtx) file and convert it to a PyTorch sparse COO tensor.

    Args:
        file_path (str): Path to the Matrix Market file.

    Returns:
        torch.Tensor: A PyTorch sparse tensor in COO format.
    """
    coo = mmread(file_path).tocoo()
    indices = torch.tensor(np.vstack((coo.row, coo.col)), dtype=torch.int64)
    values = torch.tensor(coo.data, dtype=torch.float32)

    return torch.sparse_coo_tensor(indices, values, coo.shape)


def coo_to_torch(coo, dense=True):
    """
    Convert a SciPy COO sparse matrix to a PyTorch tensor.

    Args:
        coo (scipy.sparse.coo_matrix): Input sparse matrix.
        dense (bool): If True, return a dense tensor. Default is True.

    Returns:
        torch.Tensor: Converted PyTorch tensor.
    """
    coo = np.ascontiguousarray(coo)
    indices = torch.tensor(np.vstack((coo.row, coo.col)), dtype=torch.int64)
    values = torch.tensor(coo.data, dtype=torch.float32)
    tensor = torch.sparse_coo_tensor(indices, values, coo.shape)

    return tensor.to_dense() if dense else tensor


def csc_to_torch_tensor(csc_matrix, to_dense=False):
    """
    Convert a SciPy CSC sparse matrix to a PyTorch tensor.

    Args:
        csc_matrix (scipy.sparse.csc_matrix): Input CSC sparse matrix.
        to_dense (bool): Return dense tensor if True.

    Returns:
        torch.Tensor: PyTorch tensor (dense or sparse).
    """
    coo = csc_matrix.tocoo()
    indices = torch.tensor([coo.row, coo.col], dtype=torch.int64)
    values = torch.tensor(coo.data, dtype=torch.float32)

    tensor = torch.sparse_coo_tensor(indices, values, coo.shape, dtype=torch.float64)
    return tensor.to_dense() if to_dense else tensor


def zero_pad(tensor, tile_size):
    """
    Pads the last two dimensions of a tensor with zeros so they are multiples of tile_size.

    Args:
        tensor (torch.Tensor): Input tensor of any shape.
        tile_size (int): Desired tile size.

    Returns:
        torch.Tensor: Padded tensor.
    """
    *leading_dims, height, width = tensor.shape

    pad_h = (tile_size - height % tile_size) % tile_size
    pad_w = (tile_size - width % tile_size) % tile_size

    # Split padding equally on both sides
    pad_top = pad_h // 2
    pad_bottom = pad_h - pad_top
    pad_left = pad_w // 2
    pad_right = pad_w - pad_left

    # Note: torch.nn.functional.pad pads starting from the last dimension
    padding = (pad_left, pad_right, pad_top, pad_bottom)  # (left, right, top, bottom)

    padded_tensor = F.pad(tensor, padding, mode="constant", value=0)
    return padded_tensor


def make_8conn_adjacency(height, width, device="cpu"):
    """
    Generate 8-connected spatial adjacency matrix for a 2D grid (H x W).

    Returns:
    - indices: (2, N_edges) LongTensor of edge connections (row, col)
    - values: (N_edges,) FloatTensor of edge weights (1.0)
    - shape: (N_bins, N_bins)
    """
    H, W = height, width
    N = H * W

    def idx(y, x):
        return y * W + x  # 2D to 1D index

    edges = []

    for y in range(H):
        for x in range(W):
            center = idx(y, x)

            # 8 directions (including diagonals)
            for dy in [-1, 0, 1]:
                for dx in [-1, 0, 1]:
                    ny, nx = y + dy, x + dx
                    if 0 <= ny < H and 0 <= nx < W:
                        neighbor = idx(ny, nx)
                        edges.append((center, neighbor))

    # Convert to tensor indices
    row_idx = torch.tensor([i for i, j in edges], dtype=torch.long, device=device)
    col_idx = torch.tensor([j for i, j in edges], dtype=torch.long, device=device)
    values = torch.ones(len(edges), dtype=torch.float32, device=device)

    indices = torch.stack([row_idx, col_idx])  # shape: (2, N_edges)
    shape = (N, N)

    # Create sparse adjacency matrix
    spatial_adj = torch.sparse_coo_tensor(indices, values, shape)

    return spatial_adj.coalesce()  # returns a torch.sparse.FloatTensor


def tile_matrix(matrix, tile_size):
    """
    Split a matrix into square tiles.

    Args:
        matrix (torch.Tensor): Input tensor.
        tile_size (int): Size of each tile.

    Returns:
        torch.Tensor: Tensor of shape (num_tiles, tile_size, tile_size).
    """
    tiles_h = torch.stack(matrix.split(tile_size, dim=-1), dim=2)
    tiles = torch.stack(tiles_h.split(tile_size, dim=1), dim=3)
    return tiles.reshape(-1, tile_size, tile_size)


def tile_tensor(tensor, tile_size):
    """
    Split the last two dimensions of a tensor into square tiles.

    Args:
        tensor (torch.Tensor): Input tensor of any shape.
        tile_size (int): Size of each tile.

    Returns:
        torch.Tensor: Tensor of shape (*leading_dims, num_tiles_h, num_tiles_w, tile_size, tile_size),
                      where leading_dims are the original leading dimensions.
    """
    *leading_dims, height, width = tensor.shape
    assert (
        height % tile_size == 0 and width % tile_size == 0
    ), "Tensor dimensions must be divisible by tile_size."

    num_tiles_h = height // tile_size
    num_tiles_w = width // tile_size

    # Reshape into tiles
    tensor = tensor.reshape(
        *leading_dims, num_tiles_h, tile_size, num_tiles_w, tile_size
    )

    # Correct permutation:
    # We want: (*leading_dims, num_tiles_h, num_tiles_w, tile_size, tile_size)
    dims = list(range(len(leading_dims))) + [
        len(leading_dims),
        len(leading_dims) + 2,
        len(leading_dims) + 1,
        len(leading_dims) + 3,
    ]
    tensor = tensor.permute(*dims)

    return tensor


def get_tile(matrix, index, tile_size):
    """
    Extract a tile from a matrix or batch of matrices.

    Args:
        matrix (torch.Tensor): Input matrix or batch.
        index (int): Tile index.
        tile_size (int): Size of the tile.

    Returns:
        torch.Tensor: Extracted tile.
    """
    n = matrix.size(1) if matrix.ndim == 2 else matrix.size(2)
    row_start = tile_size * (index * tile_size // n)
    col_start = (index * tile_size) % n

    if matrix.ndim > 2:
        return matrix[
            :, row_start : row_start + tile_size, col_start : col_start + tile_size
        ]
    return matrix[row_start : row_start + tile_size, col_start : col_start + tile_size]


def get_landmarks_signal(level, tile_sizes, landmark_matrices, index=None):
    """
    Compute landmark signal at a given level.

    Args:
        level (int): Current resolution level.
        tile_sizes (list[int]): Tile sizes at each level.
        landmark_matrices (torch.Tensor): Batched matrices.
        index (int, optional): Tile index if level > 0.

    Returns:
        torch.Tensor: Normalized landmark signal.
    """
    if level > 0 and index is None:
        raise ValueError("Index must be specified for level > 0.")

    tile_fn = torch.vmap(tile_matrix)
    tile_size = tile_sizes[level]

    if level == 0:
        tiles = tile_fn(landmark_matrices, tile_size=tile_size)
    else:
        tile = get_tile(landmark_matrices, index, tile_sizes[level - 1])
        tiles = tile_fn(tile, tile_size=tile_size)

    signal = tiles.sum(dim=(-1, -2))
    signal = signal / signal.sum(dim=-1, keepdim=True)
    return torch.nan_to_num(signal)


def change_landmark_resolution(tile_size, landmark_matrices):
    """
    Adjust resolution of landmark matrices using tile size.

    Args:
        tile_size (int): Desired tile size.
        landmark_matrices (torch.Tensor): Input matrices.

    Returns:
        torch.Tensor: Rescaled landmark signal.
    """
    tile_fn = torch.vmap(tile_matrix)
    tiles = tile_fn(landmark_matrices, tile_size=tile_size)
    signal = tiles.sum(dim=(-1, -2))
    return torch.nan_to_num(signal)


def normalize_likelihood(likelihood, dim=1, lambda_param=0):
    """
    Normalize likelihood using sparsemax.

    Args:
        likelihood (torch.Tensor): Input tensor.
        dim (int): Axis along which to normalize.
        lambda_param (float): Parameter for sparsemax.

    Returns:
        torch.Tensor: Normalized tensor.
    """
    output = torch.zeros_like(likelihood.T if dim == 0 else likelihood)

    for i in range(output.shape[1]):
        col = likelihood[:, i] if dim == 1 else likelihood.T[:, i]
        output[:, i] = sparsemax(col, lambda_param=lambda_param)

    output[:, -1] = 0
    return output.T if dim == 0 else output


def apply_mask(matrix, mask, return_dense=False):
    """
    Apply a mask to a sparse matrix.

    Args:
        matrix (scipy.sparse.coo_matrix): Input matrix.
        mask (scipy.sparse.coo_matrix): Mask matrix.
        return_dense (bool): Return dense matrix if True.

    Returns:
        Union[scipy.sparse.coo_matrix, np.ndarray]: Masked result.
    """
    masked_data = np.asarray(matrix[mask.row, mask.col]).reshape(-1)
    masked = scipy.sparse.coo_matrix(
        (masked_data, (mask.row, mask.col)), shape=mask.shape
    )

    return masked.todense() if return_dense else masked


def get_expectation(expression, probabilities):
    """
    Compute expected values from expression and probabilities.

    Args:
        expression (torch.Tensor): Expression matrix.
        probabilities (torch.Tensor): Probability matrix.

    Returns:
        torch.Tensor: Expected values.
    """
    return expression @ probabilities


def gaussian_filter(input_tensor, sigma=1):
    """
    Apply Gaussian filter to each channel in a tensor.

    Args:
        input_tensor (torch.Tensor): Input tensor.
        sigma (float): Gaussian sigma.

    Returns:
        torch.Tensor: Filtered tensor.
    """
    output = input_tensor.clone()

    for i in range(output.shape[0]):
        max_val = input_tensor[i].max()
        filtered = skimage.filters.gaussian(input_tensor[i], sigma=sigma)
        output[i] = filtered * (max_val / filtered.max()) if max_val > 0 else filtered

    return output


def get_threshold(input_tensor, k=1024):
    """
    Compute a threshold from input tensor using linear models.

    Args:
        input_tensor (torch.Tensor): Input values.
        k (int): Range around midpoint to fit models.

    Returns:
        float: Calculated threshold.
    """
    s = input_tensor[input_tensor > 0].unique()
    n = s.shape[0]
    m = int(np.floor(n / 2))

    x = np.arange
    model_1 = sklearn.linear_model.LinearRegression().fit(
        x(m - k, m + k)[:, None], s[m - k : m + k]
    )
    model_2 = sklearn.linear_model.LinearRegression().fit(
        x(n - k, n)[:, None], s[n - k : n]
    )

    x_val = (model_2.intercept_ - model_1.intercept_) / (model_1.coef_ - model_2.coef_)
    return s[int(x_val)]


def row_normalize(tensor):
    """
    Normalize each row to sum to 1.

    Args:
        tensor (torch.Tensor): Input tensor.

    Returns:
        torch.Tensor: Row-normalized tensor.
    """
    row_sums = tensor.sum(dim=1, keepdim=True)
    row_sums[row_sums == 0] = 1.0
    return tensor / row_sums


def sparsemax(z, lambda_param=0):
    """
    Sparsemax activation function with optional lambda adjustment.

    Args:
        z (torch.Tensor): Input tensor.
        lambda_param (float): Lambda parameter.

    Returns:
        torch.Tensor: Sparsemax output.
    """
    if torch.norm(z, p=1) <= 1:
        z = z / (torch.norm(z, p=1) - sys.float_info.min * sys.float_info.epsilon)

    z_sorted, _ = torch.sort(z, descending=True)
    cumsum_z = torch.cumsum(z_sorted, dim=0)
    k = torch.arange(1, z.size(0) + 1, dtype=z.dtype, device=z.device)
    condition = 1 + k * z_sorted > cumsum_z
    k_z = torch.max(k[condition]).long()
    tau = (cumsum_z[k_z - 1] - 1) / k_z

    return torch.clamp(z - tau, min=0)


def sparsegen(z, g=lambda x: x, lambda_param=0):
    """
    Generalized sparse activation function.

    Args:
        z (torch.Tensor): Input tensor.
        g (callable): Transformation function.
        lambda_param (float): Regularization parameter.

    Returns:
        torch.Tensor: Sparsegen output.
    """
    g_z = g(z)
    norm = torch.norm(g_z, p=1)
    if norm <= 1:
        g_z = g_z / (norm - sys.float_info.min * sys.float_info.epsilon)

    lambda_min = 1 - norm
    lambda_param = max(lambda_param, lambda_min)

    return sparsemax(g_z / (1 - lambda_param))


def _gaussian_kernel1d(
    sigma_bins: float, truncate: float = 3.0, device=None, dtype=None
):
    if sigma_bins <= 0:
        return torch.tensor([1.0], device=device, dtype=dtype)
    radius = int(truncate * sigma_bins + 0.5)
    x = torch.arange(-radius, radius + 1, device=device, dtype=dtype)
    k = torch.exp(-(x**2) / (2 * sigma_bins**2))
    k /= k.sum()
    return k


@torch.no_grad()
def smooth_probs_gaussian_torch_chunked(
    probs: torch.Tensor,  # (K, H, W), per-pixel sums to 1
    sigma_um: float,
    bin_size_um: float = 2.56,
    truncate: float = 3.0,
    mask: torch.Tensor | None = None,  # (H, W) 0/1
    eps: float = 1e-12,
    chunk_k: int = 32,
    device: torch.device | None = None,  # where to compute (cpu or cuda)
    inplace: bool = False,
    pad_mode: str = "reflect",
):
    """
    Memory-efficient Gaussian smoothing on the simplex:
    - Assumes probs of shape (K, H, W).
    - Processes channels in chunks of size `chunk_k`.
    - Two passes: (1) smooth & accumulate sum over K, (2) renormalize per pixel.
    - Denominator conv(mask) computed once and reused.
    """
    assert probs.ndim == 3, "probs must be (K, H, W)"
    K, H, W = probs.shape
    if device is None:
        device = probs.device  # compute on same device by default

    dtype = torch.float32 if probs.dtype == torch.float64 else probs.dtype
    sigma_bins = max(sigma_um / bin_size_um, 0.0)

    # Build separable kernels
    kx = _gaussian_kernel1d(sigma_bins, truncate, device=device, dtype=dtype).view(
        1, 1, 1, -1
    )
    ky = _gaussian_kernel1d(sigma_bins, truncate, device=device, dtype=dtype).view(
        1, 1, -1, 1
    )
    pad_x = (kx.shape[-1] - 1) // 2
    pad_y = (ky.shape[-2] - 1) // 2

    def conv_sep(Z: torch.Tensor):
        Z = F.pad(Z, (pad_x, pad_x, pad_y, pad_y), mode=pad_mode)
        Z = F.conv2d(Z, kx)
        Z = F.conv2d(Z, ky)
        return Z

    # Prepare mask and its convolution (denominator) ONCE
    if mask is None:
        M = torch.ones((1, 1, H, W), device=device, dtype=dtype)
    else:
        M = mask.to(device=device, dtype=dtype).clamp(0, 1).unsqueeze(0).unsqueeze(0)

    denom = conv_sep(M).clamp_min(1e-8)  # (1,1,H,W)

    # Output buffer (can write in-place to save RAM)
    out = (
        probs
        if inplace
        else torch.empty_like(probs, device=probs.device, dtype=probs.dtype)
    )

    # First pass: smooth chunks, accumulate per-pixel sum over K
    sum_map = torch.zeros((H, W), device=device, dtype=dtype)
    for c0 in range(0, K, chunk_k):
        c1 = min(K, c0 + chunk_k)

        # Move chunk to compute device
        X = probs[c0:c1].to(device=device, dtype=dtype)  # (C, H, W)
        X = X.unsqueeze(1)  # (C,1,H,W)

        num = conv_sep(X * M)  # (C,1,H,W)
        Y = (num / denom).squeeze(1)  # (C, H, W)

        # Accumulate per-pixel total (for later normalization)
        sum_map += Y.sum(dim=0)  # sum over channel axis -> (H, W)

        # Write chunk back to output (on original device)
        out[c0:c1] = Y.to(out.device, dtype=out.dtype)

        # free intermediates
        del X, num, Y
        if device.type == "cuda":
            torch.cuda.empty_cache()

    # Second pass: renormalize to the simplex per pixel
    sum_map = sum_map.clamp_min(eps).unsqueeze(0)  # (1, H, W)
    for c0 in range(0, K, chunk_k):
        c1 = min(K, c0 + chunk_k)
        out[c0:c1] = out[c0:c1] / sum_map.to(out.device, dtype=out.dtype)

    return out


def find_all_local_maxima(
    prob_map: torch.Tensor,
    threshold: float | None = None,
    one_per_plateau: bool = False,
):
    """
    Find all local maxima in a 2D tensor (8-connected neighborhood).

    Args:
        prob_map: 2D tensor (H, W)
        threshold: keep only maxima >= threshold (optional)
        one_per_plateau: if True, choose a single representative pixel per flat plateau

    Returns:
        coords: (N, 2) int tensor of (row, col) for each local maximum
        values: (N,) tensor of values at those coords
        mask: (H, W) bool tensor where True marks local maxima
    """
    assert prob_map.ndim == 2, "prob_map must be 2D (H, W)"
    x = prob_map.unsqueeze(0).unsqueeze(0)  # [1,1,H,W]

    # Max filter (3x3) -> local peaks are where x equals its pooled max
    pooled = F.max_pool2d(x, kernel_size=3, stride=1, padding=1)
    mask = x == pooled

    # Optional thresholding
    if threshold is not None:
        mask &= x >= threshold

    # Optional: reduce flat plateaus to one pixel (deterministic tiebreak)
    if one_per_plateau:
        H, W = prob_map.shape
        yy, xx = torch.meshgrid(
            torch.arange(H, device=prob_map.device),
            torch.arange(W, device=prob_map.device),
            indexing="ij",
        )
        # Tiny monotonic epsilon breaks ties deterministically by (row, col)
        eps = (yy * W + xx).to(prob_map.dtype) * 1e-12
        x_eps = x + eps
        pooled_eps = F.max_pool2d(x_eps, kernel_size=3, stride=1, padding=1)
        mask = mask & (x_eps == pooled_eps)

    mask2d = mask[0, 0]
    coords = torch.nonzero(mask2d, as_tuple=False)  # (N, 2) -> (row, col)
    values = prob_map[mask2d]  # (N,)
    return coords, values, mask2d
