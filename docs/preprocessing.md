<!-- Part of the KSTAT Tutorial. Return to [Main Tutorial](tutorial.md) -->

# 1. Preprocessing Pipeline

This section explains how to preprocess the raw reads and tissue images into formats used by the downstream model.

---

## 1.1 Project Directory Structure

Ensure the following structure exists:

```text
    kstat/
    |-- config/             # Experiment configuration files
    |-- data/               # Input CSVs and output matrices
    |-- image/              # Background tissue image and outputs
    |-- python/
    |   |-- kstat/          # Python modules
    |   |-- scripts/        # Execution scripts
    |-- docs/               # This tutorial
    |-- R/                  # R scripts for integration
```

## 1.2 Configuration File
A sample JSON configuration file (config/mouse_e18.5d.json) might look like:
```json
    {
		"project_root": "/absolute/path/to/kstat/",
		"experiment_name": "mouse_e18.5d",
		"stage": "e18.5d",
		"scale_factor": 16,
		 "reads_file": "cartana_1b_4_e18.5d_high.csv",
		 "original_background_image": "mouse_e18.5d_dapi.jpg",
		 "scaled_height": 926,
		 "scaled_width": 1039
    }
```
Save this in `config/mouse_e18.5d.json.`

## 1.3 Preprocessing Pipeline Overview

The preprocessing pipeline prepares the image, landmarks, and masks used for spatial analysis.

This step is implemented in two files:

    - The reusable module: python/kstat/preprocessing.py
    - The driver script: python/scripts/run_preprocessing.py

The pipeline performs:

    1. Load configuration parameters from JSON
    2. Load and clean landmark read data
    3. Process the background tissue image
    4. Generate a sparse landmark matrix
    5. Generate and resize a binary mask
    6. Resize and save the background image
    7. Save all outputs for downstream R analysis

---

### 1.3.1 Load Configuration
```python
    from kstat.preprocessing import load_config

    config = load_config("config/mouse_e18.5d.json")
```
---

### 1.3.2 Load and Clean Landmark Reads
```python
    import pandas as pd
    import numpy as np

    landmark_reads = pd.read_csv(f"{config['project_root']}data/{config['reads_file']}")
    landmark_reads[["X", "Y"]] = landmark_reads[["X", "Y"]].apply(np.floor).astype(int)

    counts_df = landmark_reads[["X", "Y"]].value_counts().reset_index(name="counts")
    landmark_reads = landmark_reads.merge(counts_df, on=["X", "Y"], how="left")
```
---

### 1.3.3 Process Background Image
```python
    from kstat.preprocessing import process_image
	padded_image_tensor, padding = process_image(
        f"{config['project_root']}image/{config['original_background_image']}",
        config
    )
```

This function:

    - Loads the input tissue image using PIL and converts it to a tensor
    - Pads the image so its height and width are divisible by scale_factor
    - Returns the padded image and the padding applied
	
### 1.3.4 Generate Landmark Matrix
```python
    from kstat.preprocessing import prepare_landmarks

    landmarks_matrix = prepare_landmarks(landmark_reads, padded_image_tensor, config)
```
This function:

    - Groups landmark reads by gene
    - Converts coordinates into a sparse tensor per gene
    - Tiles the image into scaled grid cells
    - Aggregates counts per tile
    - Returns a sparse gene-by-tile matrix

---

### 1.3.5 Generate Binary Mask and Resize
```python
    from kstat.preprocessing import apply_dilation
    import torchvision.transforms.functional as TF

    mask = padded_image_tensor.any(dim=0).int()
    dilated_mask = apply_dilation(mask)

    scaled_dilated_mask = (
        TF.resize(
            dilated_mask[None, :, :],
            (config["scaled_height"], config["scaled_width"])
        )
        .squeeze()
        .to_sparse()
    )
```
This step:

    - Detects tissue presence by creating a binary mask
    - Applies binary dilation to expand the region
    - Scales the mask to match the desired output resolution
    - Converts it to a sparse tensor for efficient storage

---

### 1.3.6 Resize and Save Image
```python
    from torchvision import transforms

    scaled_padded_image_tensor = TF.resize(
        padded_image_tensor, (config["scaled_height"], config["scaled_width"])
    )
    scaled_padded_image = transforms.ToPILImage()(scaled_padded_image_tensor)

This converts the padded image tensor back to a standard image format:

    - Downsamples the image to match the scaled dimensions
    - Converts it to a PIL image for saving and visualization
```
---

### 1.3.7 Save Outputs
```python
    from kstat.preprocessing import save_scaled_outputs, save_config

    save_scaled_outputs(
        config,
        landmarks_matrix,
        scaled_dilated_mask,
        scaled_padded_image
    )

    save_config(config, "config/mouse_e18.5d.json")
```
This function saves the following files:

    - data/mouse_e18.5d_landmarks_matrix.mtx
    - data/mouse_e18.5d_scaled_mask.mtx
    - image/mouse_e18.5d_background_scaled.jpg

[Next: Inference Pipeline](inference.md)
