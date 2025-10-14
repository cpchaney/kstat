import sys
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from cellpose_omni import io, models
from matplotlib.colors import Normalize
from scipy.io import mmread
from scipy.sparse import csr_matrix, diags
from skimage.exposure import equalize_adapthist
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries, watershed
from skimage.transform import rescale, resize
from tqdm import tqdm

# ========== CONFIG ==========
SCALE = 0.8
DIAMETER_UM = 50
SCALED_WIDTH = 1039
SCALED_HEIGHT = 926


# ========== MODULES ==========


def load_dapi_image(path, scale=0.8):
    dapi_hr = io.imread(path)
    dapi_eq = equalize_adapthist(dapi_hr, clip_limit=0.03)
    dapi_rescaled = rescale(dapi_eq, scale, preserve_range=True, anti_aliasing=True)
    return dapi_rescaled


def run_omnipose_segmentation(image, diameter_px):
    model = models.CellposeModel(model_type="nuclei", gpu=True)
    masks, flows, styles = model.eval(
        image,
        channels=[0, 0],
        diameter=diameter_px,
        do_3D=False,
        tile=True,
        augment=False,
        stitch_threshold=0,
    )
    return masks


def normalize_binwise_probabilities(prob_maps, cell_types, exclude="Und"):
    if isinstance(exclude, str):
        exclude = [exclude]

    valid_mask = ~cell_types.isin(exclude)
    valid_indices = np.where(valid_mask)[0]
    filtered_cell_types = cell_types.iloc[valid_indices]
    filtered_prob_maps = prob_maps[valid_indices, :]  # rows = snRNA-seq cells

    col_sums = np.array(filtered_prob_maps.sum(axis=0)).flatten()
    scaling_factors = np.ones_like(col_sums)
    nonzero_cols = col_sums > 0
    scaling_factors[nonzero_cols] = 1.0 / col_sums[nonzero_cols]

    D_inv = diags(scaling_factors)
    normalized_maps = filtered_prob_maps @ D_inv
    return normalized_maps.tocsr(), filtered_cell_types


def expand_masks_with_watershed(masks_rescaled, prob_maps, shape):
    pixel_max_probs = prob_maps.max(axis=0).toarray().flatten()
    pixel_max_probs_img = pixel_max_probs.reshape(shape)

    topography = -pixel_max_probs_img
    tissue_mask = pixel_max_probs_img > 0

    cell_masks = watershed(topography, markers=masks_rescaled, mask=tissue_mask)
    return cell_masks


def aggregate_cell_type_probabilities(cell_masks, prob_maps, cell_types):
    H, W = cell_masks.shape
    num_cell_types = prob_maps.shape[0]
    num_cells = cell_masks.max()
    cell_probs = np.zeros((num_cells, num_cell_types))

    for i, region in tqdm(enumerate(regionprops(cell_masks)), total=num_cells):
        coords = region.coords
        flat_indices = coords[:, 0] * W + coords[:, 1]
        prob_submatrix = prob_maps[:, flat_indices]
        cell_probs[i] = np.array(prob_submatrix.sum(axis=1)).flatten()
        prob_mass = cell_probs[i].sum()
        if prob_mass > 0:
            cell_probs[i] = cell_probs[i] / prob_mass

    return cell_probs


def build_df_cells(cell_probs, cell_types):
    unique_types = cell_types.unique()
    unique_types.sort()

    type_to_indices = {
        cell_type: np.where(cell_types.values == cell_type)[0]
        for cell_type in unique_types
    }

    num_cells, _ = cell_probs.shape
    num_types = len(unique_types)

    cell_type_probs = np.zeros((num_cells, num_types))

    for j, cell_type in enumerate(unique_types):
        idx = type_to_indices[cell_type]
        cell_type_probs[:, j] = cell_probs[:, idx].sum(axis=1)

    assigned_indices = np.argmax(cell_type_probs, axis=1)
    assigned_cell_types = [unique_types[i] for i in assigned_indices]

    df_cells = pd.DataFrame(
        {
            "cell_id": np.arange(1, num_cells + 1),
            "assigned_cell_type": assigned_cell_types,
        }
    )

    for j, cell_type in enumerate(unique_types):
        df_cells[f"{cell_type}_prob"] = cell_type_probs[:, j]

    return df_cells, unique_types


def plot_cell_overlay(
    dapi,
    cell_masks,
    df_cells,
    target_type,
    mode="assigned",
    cmap_name="plasma",
    save_path=None,
):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(dapi, cmap="gray")

    if mode == "assigned":
        selected_ids = df_cells[df_cells["assigned_cell_type"] == target_type][
            "cell_id"
        ]
        mask_selected = np.isin(cell_masks, selected_ids)
        overlay = np.ma.masked_where(~mask_selected, mask_selected)
        ax.imshow(overlay, cmap="Greens", alpha=0.5)
        ax.set_title(f"Cells assigned to: {target_type}")

    elif mode == "likelihood":
        prob_column = f"{target_type}_prob"
        if prob_column not in df_cells.columns:
            raise ValueError(f"{prob_column} not found in df_cells.")
        cell_id_to_prob = dict(zip(df_cells["cell_id"], df_cells[prob_column]))
        heatmap = np.zeros_like(cell_masks, dtype=np.float32)
        for cell_id, prob in cell_id_to_prob.items():
            if prob > 0:
                heatmap[cell_masks == cell_id] = prob
        norm = Normalize(vmin=0, vmax=np.max(df_cells[prob_column]))
        cmap = cm.get_cmap(cmap_name)
        colored_overlay = cmap(norm(heatmap))
        ax.imshow(colored_overlay, alpha=0.6)
        ax.set_title(f"{target_type} cell-type likelihoods")

    else:
        raise ValueError("mode must be 'assigned' or 'likelihood'")

    ax.axis("off")
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


# ========== MAIN PIPELINE ==========


def main():
    with open("../../config/config.yaml") as f:
        cfg = yaml.safe_load(f)

    PROJECT_ROOT = Path(cfg["project_root"])

    dapi_rescaled = load_dapi_image(PROJECT_ROOT / "image/mouse_e18.5d_dapi.tif", SCALE)
    diameter_px = int(DIAMETER_UM * SCALE)
    masks = run_omnipose_segmentation(dapi_rescaled, diameter_px)

    # Resize masks to match probability map shape
    masks_rescaled = resize(
        masks,
        output_shape=(SCALED_HEIGHT, SCALED_WIDTH),
        order=0,
        preserve_range=True,
        anti_aliasing=False,
    ).astype(np.int32)

    # Load probability maps and cell type table
    prob_maps = mmread(
        PROJECT_ROOT / "data/mouse_e18.5d_cells_bins_probabilities.mtx"
    ).tocsr()

    cell_type_table = pd.read_csv(
        PROJECT_ROOT / "../kidney_development/mouse_e18.5d/data/cell_type_table.csv",
        names=["cell_id", "coarse_cell_type", "cell_type"],
        skiprows=1,
    )
    normalized_prob_maps, filtered_cell_types = normalize_binwise_probabilities(
        prob_maps, cell_type_table["cell_type"], exclude="Und"
    )
    unique_bins_key = pd.read_csv(
        PROJECT_ROOT / "data/mouse_e18.5d_unique_bins_key.csv"
    )
    bin_indices = unique_bins_key["index"].to_numpy(dtype=np.int64)

    cell_masks = expand_masks_with_watershed(
        masks_rescaled,
        normalized_prob_maps[:, bin_indices],
        (SCALED_HEIGHT, SCALED_WIDTH),
    )
    cell_probs = aggregate_cell_type_probabilities(
        cell_masks, normalized_prob_maps[:, bin_indices], filtered_cell_types
    )
    df_cells, unique_types = build_df_cells(cell_probs, filtered_cell_types)

    # Example: show visualizations for one target type
    target_type = "Int_8"
    plot_cell_overlay(dapi_rescaled, cell_masks, df_cells, target_type, mode="assigned")
    plot_cell_overlay(
        dapi_rescaled, cell_masks, df_cells, target_type, mode="likelihood"
    )


if __name__ == "__main__":
    main()
