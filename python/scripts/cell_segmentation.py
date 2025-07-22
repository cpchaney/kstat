import sys
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from cellpose_omni import io, models, utils
from matplotlib import colors
from matplotlib.colors import Normalize
from scipy.io import mmread
from scipy.ndimage import distance_transform_edt
from scipy.sparse import csr_matrix, diags
from skimage import exposure
from skimage.color import rgb2gray
from skimage.exposure import equalize_adapthist
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries, watershed
from skimage.transform import rescale, resize
from tqdm import tqdm


def normalize_binwise_probabilities(prob_maps, cell_types, exclude="Und"):
    """
    Remove excluded cell types and normalize columns (spatial bins) so they sum to 1.
    Keeps zero columns zero.

    Parameters:
        prob_maps (csr_matrix): Shape (num_cells, num_bins), probability per cell/bin
        cell_types (pd.Series): Length = num_cells (rows of prob_maps)
        exclude (str or list): Cell type(s) to exclude, e.g., "Und"

    Returns:
        normalized_maps (csr_matrix): Shape (filtered_cells, num_bins), normalized per column
        filtered_cell_types (pd.Series): Cell types used (no excluded types)
    """
    # Convert exclude to list
    if isinstance(exclude, str):
        exclude = [exclude]

    # Step 1: filter out excluded cell types (rows)
    valid_mask = ~cell_types.isin(exclude)
    valid_indices = np.where(valid_mask)[0]
    filtered_cell_types = cell_types.iloc[valid_indices]
    filtered_prob_maps = prob_maps[
        valid_indices, :
    ]  # shape: (filtered_cells, num_bins)

    # Step 2: normalize each column (bin)
    col_sums = np.array(filtered_prob_maps.sum(axis=0)).flatten()  # shape: (num_bins,)
    nonzero_cols = col_sums > 0
    scaling_factors = np.ones_like(col_sums)
    scaling_factors[nonzero_cols] = 1.0 / col_sums[nonzero_cols]

    # Build diagonal matrix to apply scaling
    D_inv = diags(scaling_factors)  # shape: (num_bins, num_bins)

    # Apply normalization: each column sums to 1 (or stays 0)
    normalized_maps = filtered_prob_maps @ D_inv

    return normalized_maps.tocsr(), filtered_cell_types


def plot_assigned_cells(dapi, cell_masks, df_cells, target_type, color="lime"):
    mask_selected = np.isin(
        cell_masks, df_cells[df_cells["assigned_cell_type"] == target_type]["cell_id"]
    )

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(dapi, cmap="gray")

    # Overlay mask as a binary outline or filled
    overlay = np.ma.masked_where(~mask_selected, mask_selected)
    ax.imshow(overlay, cmap=plt.cm.get_cmap("Greens"), alpha=0.5)

    ax.set_title(f"Cells assigned to: {target_type}")
    ax.axis("off")
    plt.tight_layout()
    plt.show()


def plot_cell_type_likelihoods(
    dapi, cell_masks, df_cells, target_type, cmap_name="plasma"
):
    prob_column = f"{target_type}_prob"

    if prob_column not in df_cells.columns:
        raise ValueError(f"{prob_column} not found in df_cells.")

    # Map cell_id ? probability for the target type
    cell_id_to_prob = dict(zip(df_cells["cell_id"], df_cells[prob_column]))

    # Create an image with float values for each pixel in cell_masks
    heatmap = np.zeros_like(cell_masks, dtype=np.float32)

    for cell_id, prob in cell_id_to_prob.items():
        if prob > 0:
            heatmap[cell_masks == cell_id] = prob

    # Normalize to [0, 1] for colormap mapping
    norm = Normalize(vmin=0, vmax=np.max(df_cells[prob_column]))
    cmap = cm.get_cmap(cmap_name)

    colored_overlay = cmap(norm(heatmap))

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(dapi, cmap="gray")
    ax.imshow(colored_overlay, alpha=0.6)
    ax.set_title(f"{target_type} cell-type likelihoods")
    ax.axis("off")
    plt.tight_layout()
    plt.show()


# Use current working directory as fallback
script_dir = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()

with open("../../config/config.yaml") as f:
    cfg = yaml.safe_load(f)

PROJECT_ROOT = Path(cfg["project_root"])

# Load your DAPI image
dapi_hr = io.imread(
    PROJECT_ROOT / "image/mouse_e18.5d_dapi.tif"
)  # shape: (H, W) or (N, H, W)

dapi_eq = equalize_adapthist(dapi_hr, clip_limit=0.03)

# Downscale by 0.8x (or more if needed)
scale = 0.8
dapi_rescaled = rescale(dapi_eq, scale, preserve_range=True, anti_aliasing=True)

# Adjust diameter accordingly
scaled_diameter = int(50 * scale)

# Load the Omnipose model
model = models.CellposeModel(model_type="nuclei", gpu=True)

# Run segmentation
masks, flows, styles = model.eval(
    dapi_rescaled,
    channels=[0, 0],
    diameter=scaled_diameter,
    do_3D=False,
    tile=True,
    augment=False,
    stitch_threshold=0,
)

cell_type_table = pd.read_csv(
    PROJECT_ROOT / "../kidney_development/mouse_e18.5d/data/cell_type_table.csv",
    names=["cell_id", "coarse_cell_type", "cell_type"],
    skiprows=1,
)

cell_types = cell_type_table["cell_type"]

# Remove "Und" cells from both prob_maps and cell_types
valid_indices = np.where(cell_types.values != "Und")[0]
filtered_cell_types = cell_types.iloc[valid_indices]

unique_bins_key = pd.read_csv(PROJECT_ROOT / "data/mouse_e18.5d_unique_bins_key.csv")
prob_maps = mmread(PROJECT_ROOT / "data/mouse_e18.5d_cells_bins_probabilities.mtx")
prob_maps = prob_maps.tocsr()

normalized_prob_maps, filtered_cell_types = normalize_binwise_probabilities(
    prob_maps, cell_type_table["cell_type"], exclude="Und"
)

normalized_prob_maps = normalized_prob_maps[:, unique_bins_key["index"]]

SCALED_WIDTH = 1039
SCALED_HEIGHT = 926

# Suppose you downscaled by factor S (e.g., 0.1 or similar)
scale_y = SCALED_HEIGHT / dapi_rescaled.shape[1]
scale_x = SCALED_WIDTH / dapi_rescaled.shape[0]

# Resize masks using nearest-neighbor to preserve labels
masks_rescaled = resize(
    masks,
    output_shape=(SCALED_HEIGHT, SCALED_WIDTH),  # (H, W)
    order=0,  # nearest-neighbor
    preserve_range=True,
    anti_aliasing=False,
).astype(np.int32)

########################################
# Create Per-Pixel Max-Probability Map #
########################################
# Create per-pixel "confidence" map (highest probability across cell types)
pixel_max_probs = (
    normalized_prob_maps.max(axis=0).toarray().flatten()
)  # shape: (H * W,)
pixel_max_probs_img = pixel_max_probs.reshape((SCALED_HEIGHT, SCALED_WIDTH))

# Allow expansion only where there is probability signal
tissue_mask = pixel_max_probs_img > 0  # or use a binary tissue mask

# Invert for watershed: high probability ? low cost
topography = -pixel_max_probs_img

cell_masks = watershed(topography, markers=masks_rescaled, mask=tissue_mask)

H, W = cell_masks.shape
num_cell_types = normalized_prob_maps.shape[0]
num_cells = cell_masks.max()
cell_probs = np.zeros((num_cells, num_cell_types))

for i, region in tqdm(enumerate(regionprops(cell_masks))):
    coords = region.coords
    rows, cols = coords[:, 0], coords[:, 1]
    flat_indices = rows * W + cols

    prob_submatrix = normalized_prob_maps[:, flat_indices]
    cell_probs[i] = np.array(prob_submatrix.mean(axis=1)).flatten()

unique_types = filtered_cell_types.unique()
unique_types.sort()

# Map cell types to indices
type_to_indices = {
    cell_type: np.where(filtered_cell_types.values == cell_type)[0]
    for cell_type in unique_types
}

num_cells, _ = cell_probs.shape
num_types = len(unique_types)

# New array: each row is per segmented cell, each column is a cell type
cell_type_probs = np.zeros((num_cells, num_types))

# Aggregate by mean or sum over snRNA-seq cells of each type
for j, cell_type in enumerate(unique_types):
    idx = type_to_indices[cell_type]
    cell_type_probs[:, j] = cell_probs[:, idx].sum(axis=1)  # or use .sum(axis=1)

assigned_type_indices = np.argmax(cell_type_probs, axis=1)
assigned_cell_types = [unique_types[i] for i in assigned_type_indices]

# Get unique values and their counts
unique_values, counts = np.unique(assigned_cell_types, return_counts=True)

# Convert to DataFrame
counts_df = pd.DataFrame({"value": unique_values, "count": counts})

print(counts_df)

df_cells = pd.DataFrame(
    {
        "cell_id": np.arange(1, num_cells + 1),
        "assigned_cell_type": assigned_cell_types,
    }
)

# Add full probability vector
for j, cell_type in enumerate(unique_types):
    df_cells[f"{cell_type}_prob"] = cell_type_probs[:, j]


#################
# Visualization #
#################
# Nuclei segmentation
# If the image is RGB, use it directly.
# If it's grayscale (shape: H x W), you can use it with cmap='gray'
fig, ax = plt.subplots(figsize=(10, 10))

# Show the base image
if dapi_hr.ndim == 3:
    ax.imshow(dapi_rescaled)  # RGB image
else:
    ax.imshow(dapi_rescaled, cmap="gray")  # Grayscale image

# Overlay the segmentation masks
masked = np.ma.masked_where(masks == 0, masks)

# Random colormap to show different cells
cmap = plt.cm.nipy_spectral.copy()
cmap.set_bad(color="none")  # So 0s (background) are transparent
ax.imshow(masked, cmap=cmap, alpha=0.5)

# ax.set_title("Omnipose Cell Segmentation Overlay")
# ax.axis('off')
# plt.tight_layout()
# plt.show()

ax.axis("off")
plt.tight_layout()
fig.savefig(
    PROJECT_ROOT / "image/cell_segmentation_overlay.png", bbox_inches="tight", dpi=300
)
plt.close(fig)

# Watershed masks overlay
fig, ax = plt.subplots(figsize=(10, 10))

# Show DAPI image
ax.imshow(dapi_rescaled, cmap="gray")

# Mask out background (label 0) so we can overlay cells
masked_cells = np.ma.masked_where(cell_masks == 0, cell_masks)

# Use a colorful colormap to show cell territories
cmap = plt.cm.nipy_spectral.copy()
cmap.set_bad(color="none")  # makes background transparent

# Overlay the segmented cells
ax.imshow(masked_cells, cmap=cmap, alpha=0.5)

ax.axis("off")
plt.tight_layout()
fig.savefig(
    PROJECT_ROOT / "image/cell_segmentation_watershed_masks_overlay.png",
    bbox_inches="tight",
    dpi=300,
)
plt.close(fig)

# Watershed mask boundaries
boundaries = find_boundaries(cell_masks, mode="outer")

fig, ax = plt.subplots(figsize=(10, 10))

# Base DAPI image
ax.imshow(dapi_rescaled, cmap="gray")

# Overlay boundaries
ax.imshow(boundaries, cmap="hot", alpha=0.8)

# Clean up axes and layout
ax.axis("off")
plt.tight_layout()

# Save to file
fig.savefig(
    PROJECT_ROOT / "image/cell_segmentation_watershed_masks_boundaries.png",
    bbox_inches="tight",
    dpi=300,
)

# Close the figure to free memory
plt.close(fig)

# Assigned cell types to segmented cells
target_type = "Int_9"

# Show assigned cells
plot_assigned_cells(dapi_rescaled, cell_masks, df_cells, target_type)

# Show likelihoods
plot_cell_type_likelihoods(dapi_rescaled, cell_masks, df_cells, target_type)
