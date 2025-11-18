#!/usr/bin/env python
"""
Preprocessing pipeline for spatial landmark data and background image.

Usage:
    python scripts/run_preprocessing.py config/mouse_e18.5d.json
"""

import sys

import numpy as np
import pandas as pd
import torchvision.transforms as transforms
import torchvision.transforms.functional as TF
from kstat.preprocessing import (apply_dilation, load_config,
                                 prepare_landmarks, process_image, save_config,
                                 save_scaled_outputs)


def main(config_path: str):
    # ---------------------------------------------------
    # Step 1: Load configuration
    # ---------------------------------------------------
    config = load_config(config_path)

    # ---------------------------------------------------
    # Step 2: Load landmark reads and count per pixel
    # ---------------------------------------------------
    landmark_reads = pd.read_csv(f"{config['project_root']}data/{config['reads_file']}")

    # Round (X, Y) to integer pixel coordinates
    landmark_reads[["X", "Y"]] = landmark_reads[["X", "Y"]].apply(np.floor).astype(int)

    # Count number of reads per pixel
    counts_df = landmark_reads[["X", "Y"]].value_counts().reset_index(name="counts")

    # Merge counts back into dataframe
    landmark_reads = landmark_reads.merge(counts_df, on=["X", "Y"], how="left")

    # ---------------------------------------------------
    # Step 3: Load background image and apply padding
    # ---------------------------------------------------
    padded_image_tensor, padding = process_image(
        f"{config['project_root']}image/{config['original_background_image']}",
        config,
    )

    # ---------------------------------------------------
    # Step 4: Generate landmark matrix (pixels × genes)
    # ---------------------------------------------------
    landmarks_matrix = prepare_landmarks(landmark_reads, padded_image_tensor, config)

    # ---------------------------------------------------
    # Step 5: Generate binary mask and dilated mask
    # ---------------------------------------------------

    # Create a binary mask from padded image (1 if any channel is nonzero)
    mask = padded_image_tensor.any(dim=0).int()

    # Apply morphological dilation to the mask
    dilated_mask = apply_dilation(mask)

    # Resize dilated mask to scaled resolution and convert to sparse tensor
    scaled_dilated_mask = (
        TF.resize(
            dilated_mask[None, :, :],  # add batch channel
            (config["scaled_height"], config["scaled_width"]),
        )
        .squeeze()
        .to_sparse()
    )

    # ---------------------------------------------------
    # Step 6: Scale image and save outputs
    # ---------------------------------------------------

    # Resize padded image to scaled resolution
    scaled_padded_image_tensor = TF.resize(
        padded_image_tensor, (config["scaled_height"], config["scaled_width"])
    )

    # Convert to PIL image for saving
    scaled_padded_image = transforms.ToPILImage()(scaled_padded_image_tensor)

    # Save processed outputs (matrix, mask, image)
    save_scaled_outputs(
        config,
        landmarks_matrix,
        scaled_dilated_mask,
        scaled_padded_image,
    )

    # Save potentially updated config file
    save_config(config, config_path)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/run_preprocessing.py <path/to/config.json>")
    else:
        main(sys.argv[1])
