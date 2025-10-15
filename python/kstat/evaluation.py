import pandas as pd
import scipy.io
import torch
import tqdm
from geomloss import SamplesLoss

from kstat import model
from kstat.utils import row_normalize


def run_loocv_evaluation(
    ad,
    common_landmarks,
    landmark_counts,
    mask,
    expression_tensor,
    config,
    LAMBDA_PARAM,
    support,
    support_loc,
    spatial_expression,
):
    """
    Perform leave-one-out cross-validation (LOOCV) over landmark genes.

    For each landmark, remove it from the training set, reconstruct it
    using the remaining landmarks, and compute similarity metrics between
    predicted and observed spatial signals.

    Args:
        ad: AnnData object containing imputed expression.
        common_landmarks: List of landmark gene names.
        landmark_counts: Sparse matrix of measured landmark counts.
        mask: Binary mask for valid regions.
        expression_tensor: Tensor of imputed expression.
        config: Configuration dictionary.
        LAMBDA_PARAM: Regularization strength.
        support: Binary support vector indicating active spatial bins.
        support_loc: Coordinates of support bins.
        spatial_expression: Observed spatial landmark matrix.

    Returns:
        DataFrame with cosine similarity and Sinkhorn divergence for each gene.
    """
    loss_fn = SamplesLoss(loss="sinkhorn", p=2, blur=0.1)
    loo_similarity = []
    loo_sinkhorn = []

    for i, landmark in enumerate(tqdm.tqdm(common_landmarks, desc="LOOCV")):
        # Exclude current landmark
        loo_landmarks = common_landmarks[:i] + common_landmarks[i + 1:]

        # Get LOOCV imputed tensor
        loo_expression_tensor = (
            torch.tensor(ad[:, loo_landmarks].layers["imputed"].T)
            .float()
            .to(expression_tensor.device)
        )

        # Filter landmark counts for remaining landmarks and apply mask
        loo_masked = landmark_counts.to_dense()[
            [gene in loo_landmarks for gene in config["landmarks"]], :
        ] * mask.to_dense().reshape(-1)

        # Identify unique bin patterns and indexing key
        loo_unique_bins, loo_unique_key = torch.unique(
            (loo_masked > 0).long(), dim=1, return_inverse=True
        )
        loo_unique_bins = loo_unique_bins.to(expression_tensor.device)

        # Infer distribution parameters from expression
        loo_mean, loo_cov = model.infer_parameters_from_mixture(
            loo_expression_tensor, k=None
        )

        # Calculate bin probabilities using inferred parameters
        loo_bin_probs = model.calculate_bin_probabilities(
            loo_expression_tensor,
            loo_mean,
            loo_cov,
            loo_unique_bins,
            LAMBDA_PARAM,
        )

        # Reconstruct left-out gene
        landmark_idx = i
        recon = expression_tensor[landmark_idx] @ loo_bin_probs
        recon = recon[loo_unique_key[support.cpu()]]

        # Threshold top-k values
        k = min(
            (landmark_counts.to_dense()[landmark_idx, :] > 0).sum().item(),
            recon.shape[-1],
        )
        threshold = torch.topk(recon, k).values[-1]
        recon[recon < threshold] = 0
        recon = row_normalize(recon[None, :])

        # Compare reconstruction to ground truth
        observed = spatial_expression[landmark_idx, :]
        cos_sim = torch.nn.functional.cosine_similarity(recon, observed).item()
        sinkhorn = loss_fn(
            support_loc[recon.squeeze() > 0],
            support_loc[observed > 0],
        ).item()

        loo_similarity.append(cos_sim)
        loo_sinkhorn.append(sinkhorn)

    return pd.DataFrame(
        {
            "cosine_similarity": loo_similarity,
            "sinkhorn_divergence": loo_sinkhorn,
        },
        index=common_landmarks,
    )


def save_loocv_outputs(
    ad,
    config,
    error_landmarks,
    landmark_counts,
    mask,
    expression_tensor,
    support_loc,
    spatial_expression,
    cells_bins_probabilities,
    unique_bins_key,
):
    """
    Save LOOCV-inferred bin probability matrices for error landmarks,
    and full model reconstruction losses for all landmarks.

    Args:
        ad: AnnData object with imputed expression.
        config: Configuration dict.
        error_landmarks: Subset of landmarks with high reconstruction error.
        landmark_counts: Sparse matrix of measured landmark counts.
        mask: Binary mask of active spatial bins.
        expression_tensor: Full imputed expression tensor.
        support_loc: Spatial coordinates of support bins.
        spatial_expression: Measured spatial landmark signal.
        cells_bins_probabilities: Original model output.
        unique_bins_key: Index mapping of spatial bins.
    """
    loss_fn = SamplesLoss(loss="sinkhorn", p=2, blur=0.1)

    for landmark in tqdm.tqdm(error_landmarks, desc="Saving LOOCV matrices"):
        idx = config["common_landmarks"].index(landmark)
        loo_landmarks = (
            config["common_landmarks"][:idx] + config["common_landmarks"][idx + 1 :]
        )

        loo_expression_tensor = (
            torch.tensor(ad[:, loo_landmarks].layers["imputed"].T)
            .float()
            .to(expression_tensor.device)
        )

        loo_masked = landmark_counts.to_dense()[
            [g in loo_landmarks for g in config["landmarks"]], :
        ] * mask.to_dense().reshape(-1)

        loo_unique_bins, loo_unique_key = torch.unique(
            (loo_masked > 0).long(), dim=1, return_inverse=True
        )
        loo_unique_bins = loo_unique_bins.to(expression_tensor.device)

        loo_mean, loo_cov = model.infer_parameters_from_mixture(
            loo_expression_tensor, k=None
        )

        loo_probs = model.calculate_bin_probabilities(
            loo_expression_tensor,
            loo_mean,
            loo_cov,
            loo_unique_bins,
            config["lambda"],
        )

        # Save outputs
        base_path = f"{config['project_root']}data/{config['experiment_name']}_{landmark}_loocv"
        scipy.io.mmwrite(f"{base_path}_cells_bins_probabilities.mtx", loo_probs.cpu().to_dense())
        pd.DataFrame(loo_unique_key.cpu().numpy(), columns=["index"]).to_csv(
            f"{base_path}_unique_bins_key.csv", index=False
        )

    # Evaluate full model reconstructions
    landmark_losses = []
    projected = expression_tensor.cpu() @ cells_bins_probabilities
    projected = projected[:, unique_bins_key][..., support_loc[:, 0].cpu()]
    projected = projected.to(expression_tensor.device)

    for i, landmark in enumerate(config["common_landmarks"]):
        recon = projected[i, :]

        # Threshold top-k values
        k = min(
            (landmark_counts.to_dense()[i, :] > 0).sum().item(),
            recon.shape[-1],
        )
        threshold = torch.topk(recon, k).values[-1]
        recon[recon < threshold] = 0
        recon = row_normalize(recon[None, :])

        # Compute Sinkhorn loss
        loss = loss_fn(
            recon.reshape(-1),
            support_loc,
            spatial_expression[i, :],
            support_loc,
        ).item()
        landmark_losses.append(loss)

    # Save loss table
    pd.DataFrame(
        {"sinkhorn_divergence": landmark_losses},
        index=config["common_landmarks"],
    ).to_csv(
        f"{config['project_root']}data/{config['experiment_name']}_reconstruction_losses.csv"
    )
