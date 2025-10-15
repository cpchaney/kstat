<!-- Part of the KSTAT Tutorial. Return to [Main Tutorial](tutorial.md) -->

# 2. Model Inference and Evaluation
## 2.1 Inferring Bin-Wise Cell Probability Distributions

This section describes how to compute the probability that each reference profile corresponds to each bin in the spatial map.

The script used is:

    python/scripts/run_inference.py

This relies on:

    - kstat.model          # Gaussian mixture inference and OT
    - kstat.utils          # Normalization and I/O helpers
    - kstat.probability    # Input preparation and filtering

---

### 2.1.1 Configuration and Device
```python
    CONFIG_PATH = "config/mouse_e18.5d.json"
    device = "cuda" if torch.cuda.is_available() else "cpu"
```
This determines paths and selects GPU if available.

---

### 2.1.2 Load Landmark Data and Mask
```python
    from kstat.probability import prepare_inputs

    landmark_counts, mask, ad, config = prepare_inputs(config, device)
```
This function:

    - Loads the landmark matrix and mask
    - Loads the reference snRNA-seq dataset
    - Filters out genes with 0 expression
    - Aligns landmark genes with reference profiles

---

### 2.1.3 Prepare Bin Inputs
```python
    from kstat.probability import compute_bin_inputs

    unique_bins, support, support_loc, spatial_expression = compute_bin_inputs(
        landmark_counts, mask, landmarks, config["common_landmarks"], ad, device
    )
```
This:

    - Binarizes landmark gene activity
    - Identifies unique landmark patterns (bins)
    - Normalizes the per-bin landmark signal
    - Extracts spatial support locations

---

### 2.1.4 Infer Cell-Level Parameters
```python
    expression_tensor = torch.tensor(ad.layers["imputed"].T).float().to(device)

    mean_values, covariances = model.infer_parameters_from_mixture(
        expression_tensor, k=None
    )
```
This fits a multivariate Gaussian for each gene using the mixture model defined in `kstat.model`.

---

### 2.1.5 Compute Bin Probabilities
```python
    cells_bins_probabilities = model.calculate_bin_probabilities(
        expression_tensor,
        mean_values,
        covariances,
        unique_bins,
        LAMBDA_PARAM,
    )
```
This calculates the likelihood of a cell occupying a bin, for each cell and each bin, according to how well each reference cell matches each bin based on landmark expression distributions.

    - Output: tensor of shape (num_cells x num_bins)
    - Can be used for mapping cells to spatial locations

---

> Note: The probability matrix is large; consider saving it to disk or using only top matches.

---

### 2.1.6 Performance

    Execution time: ~X.XX seconds (varies by GPU and dataset size)

The entire pipeline runs efficiently on GPU, especially when batching is used in model calls.

## 2.2. Evaluating Reconstruction Accuracy

After inferring bin-wise probabilities, we compare the reconstructed spatial expression with the observed spatial landmarks.

The evaluation script is:

    python/scripts/run_evaluation.py

It uses:

    - Sinkhorn divergence (from Optimal Transport)
    - Cosine similarity (angular similarity between vectors)

---

### 2.2.1 Compute Expected Expression

We multiply the bin probabilities by reference expression profiles to estimate spatial expression:
```python
    expected_expression = expression_tensor @ cells_bins_probabilities
```
We normalize and subset this to only spatial support locations (where tissue is present).

---

### 2.2.2 Compare with Observed Spatial Expression

We compare the reconstructed tensor against the observed (binarized) spatial data:

    - Row-normalize both tensors
    - Match indices based on support region
    - Apply two metrics

---

### 2.2.3 Metrics Used

**Sinkhorn divergence** (from the `geomloss` package):

    - Measures alignment between two probability distributions
    - Suitable for spatial contexts (respects geometry)
    - Lower values indicate better reconstruction

**Cosine similarity**:

    - Measures angle between two vectors
    - 1.0 = perfect alignment, 0.0 = orthogonal

---

### 2.2.4 Example Output

    Mean Sinkhorn divergence: 0.0053
    Mean cosine similarity: 0.92

These metrics can help assess model quality and guide parameter tuning.

## 2.3. Leave-One-Out Cross-Validation (LOOCV)

To assess how well each landmark gene's spatial pattern can be predicted from all others, we perform leave-one-out cross-validation (LOOCV).

---

### 2.3.1 How it Works

For each landmark gene:

    - Remove the gene from the dataset
    - Train the model on all other landmarks
    - Predict spatial expression of the held-out gene
    - Compare prediction to the observed spatial distribution

This process is repeated for all landmark genes.

---

### 2.3.2 Running the Script
```bash
    python scripts/run_loocv.py config/mouse_e18.5d.json
```
This script outputs:

    - data/<experiment>_loocv_results.csv  ? LOOCV Sinkhorn + cosine per gene
    - data/<experiment>_<gene>_loocv_cells_bins_probabilities.mtx
    - data/<experiment>_<gene>_loocv_unique_bins_key.csv
    - data/<experiment>_reconstruction_losses.csv

---

### 2.3.3 Metrics

Each gene is evaluated with:

    - **Cosine similarity** between predicted and observed patterns
    - **Sinkhorn divergence** between predicted and observed distributions

Lower Sinkhorn and higher cosine indicate better predictions.

---

### 2.3.4 Identifying High-Error Genes

The top 10 genes with the highest Sinkhorn divergence are saved and their predictions exported for inspection.

This is useful for:

    - Quality control
    - Visualization
    - Identifying genes poorly modeled by current parameters
