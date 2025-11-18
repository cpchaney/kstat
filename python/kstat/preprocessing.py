import json

import numpy as np
import pandas as pd
import scipy
import torch
import torch.nn.functional as F
import tqdm
from PIL import Image
from scipy.io import mmwrite
from scipy.ndimage import binary_dilation
from torchvision import transforms


def load_config(config_path: str) -> dict:
    """Load configuration JSON file."""
    with open(config_path, "r") as file:
        return json.load(file)


def save_config(config: dict, config_path: str):
    """Save updated configuration to a JSON file."""
    with open(config_path, "w") as file:
        json.dump(config, file, indent=4)


def split_integer(n: int) -> tuple:
    """Split an integer into two nearly equal halves."""
    return (n // 2, n // 2) if n % 2 == 0 else ((n + 1) // 2, n // 2)


def calculate_padding(height: int, width: int, scale_factor: int) -> tuple:
    """Calculate padding dimensions to make dimensions divisible by scale factor."""
    vertical_padding = scale_factor - height % scale_factor
    horizontal_padding = scale_factor - width % scale_factor
    return split_integer(horizontal_padding) + split_integer(vertical_padding)


def sparse_tile_sum(sparse_tensor: torch.Tensor, tile_size: tuple) -> torch.Tensor:
    """Sum elements in sparse tensor tiles."""
    assert sparse_tensor.is_sparse, "Input tensor must be sparse"
    indices = sparse_tensor._indices()
    values = sparse_tensor._values()
    height, width = sparse_tensor.shape
    new_height = height // tile_size[0]
    new_width = width // tile_size[1]

    result = torch.zeros(new_height, new_width, dtype=values.dtype)

    for i in range(values.size(0)):
        row, col = indices[:, i]
        tile_row = row // tile_size[0]
        tile_col = col // tile_size[1]
        result[tile_row, tile_col] += values[i]

    return result


def torch_mmwrite(file_path: str, sparse_tensor: torch.Tensor):
    """Write a PyTorch sparse tensor to a Matrix Market file."""
    coo = sparse_tensor.coalesce()
    row, col = coo.indices().numpy()
    data = coo.values().numpy()
    scipy_sparse_matrix = scipy.sparse.coo_matrix((data, (row, col)), shape=coo.shape)
    mmwrite(file_path, scipy_sparse_matrix)


def prepare_landmarks(
    landmark_reads: pd.DataFrame, image_tensor: torch.Tensor, config: dict
) -> torch.Tensor:
    """Process landmark reads and generate the landmarks sparse matrix."""
    [height, width] = image_tensor.shape[1:]
    padding = calculate_padding(height, width, config["scale_factor"])
    padded_landmark_reads = landmark_reads.copy()
    padded_landmark_reads[["X", "Y"]] += [padding[0], padding[2]]

    size = (height + padding[0] + padding[1], width + padding[2] + padding[3])
    tile_size = (config["scale_factor"], config["scale_factor"])

    landmarks = []
    indices = []
    values = []

    for row_index, (gene, subframe) in enumerate(
        tqdm.tqdm(padded_landmark_reads.groupby("gene"))
    ):
        sub_indices = torch.tensor([subframe["Y"].to_numpy(), subframe["X"].to_numpy()])
        sub_values = torch.tensor(subframe["counts"].to_numpy(), dtype=torch.float32)
        sparse_matrix = torch.sparse_coo_tensor(sub_indices, sub_values, size)
        sparse_vector = (
            sparse_tile_sum(sparse_matrix, tile_size).reshape(-1).to_sparse()
        )

        landmarks.append(gene)
        indices.append(
            torch.vstack(
                (
                    torch.full((1, sparse_vector.indices().shape[1]), row_index),
                    sparse_vector.indices(),
                )
            )
        )
        values.append(sparse_vector.values())

    return torch.sparse_coo_tensor(
        torch.cat(indices, dim=1),
        torch.cat(values),
        size=(
            len(landmarks),
            (size[0] // config["scale_factor"]) * (size[1] // config["scale_factor"]),
        ),
    )


def process_image(image_path: str, config: dict) -> torch.Tensor:
    """Load and preprocess an image."""
    # Increase the allowable maximum image size (e.g., to 300 million pixels)
    Image.MAX_IMAGE_PIXELS = None  # disables the check entirely
    image = Image.open(image_path)
    transform = transforms.ToTensor()
    image_tensor = transform(image)
    [height, width] = image_tensor.shape[1:]
    padding = calculate_padding(height, width, config["scale_factor"])
    return F.pad(image_tensor, padding, "constant", 0), padding


def apply_dilation(mask: torch.Tensor, kernel_size: int = 3) -> torch.Tensor:
    """Apply binary dilation to a mask."""
    structure = np.ones((kernel_size, kernel_size), dtype=bool)
    dilated_np = binary_dilation(mask.numpy(), structure=structure)
    return torch.from_numpy(dilated_np.astype(int))


def save_scaled_outputs(
    config: dict,
    landmarks_matrix: torch.Tensor,
    scaled_dilated_mask: torch.Tensor,
    scaled_padded_image: Image.Image,
):
    """Save processed outputs."""
    torch_mmwrite(
        f"{config['project_root']}data/{config['experiment_name']}_landmarks_matrix.mtx",
        landmarks_matrix,
    )
    torch_mmwrite(
        f"{config['project_root']}data/{config['experiment_name']}_scaled_mask.mtx",
        scaled_dilated_mask,
    )
    scaled_padded_image.save(
        f"{config['project_root']}image/{config['experiment_name']}_background_scaled.jpg",
        format="JPEG",
    )
