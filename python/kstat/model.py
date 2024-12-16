import h5py
import numpy as np
import kstat.utils
import torch
import tqdm
from joblib import Parallel, delayed
from sklearn.mixture import GaussianMixture


def is_positive_semidefinite(matrix):
    try:
        np.linalg.cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False


def fit_gmm_with_best_bic(data, min_components=2, max_components=10):
    import logging

    logging.basicConfig(level=logging.INFO)

    best_bic = float("inf")
    best_model = None
    bic_scores = []

    for k in tqdm.tqdm(range(min_components, max_components + 1), desc="Fitting GMMs"):
        gmm = GaussianMixture(n_components=k, random_state=0)
        gmm.fit(data)
        bic = gmm.bic(data)
        bic_scores.append(bic)

        logging.info(f"Components: {k}, BIC: {bic}")
        if bic < best_bic:
            best_bic = bic
            best_model = gmm

    return best_model, best_bic, bic_scores


def infer_landmark_parameters(expression_tensor, index, k):
    x = expression_tensor[index, :].cpu().numpy()
    gmm = (
        GaussianMixture(n_components=k, random_state=42)
        if k
        else fit_gmm_with_best_bic(x.reshape(-1, 1), max_components=8)[0]
    )
    components = gmm.predict(x.reshape(-1, 1))
    on_component = np.argmax(gmm.means_)

    mean = torch.tensor(
        np.insert(gmm.means_[on_component], 0, x[components != on_component].mean())
    )
    covariance = torch.tensor(
        np.insert(
            gmm.covariances_[on_component], 0, x[components != on_component].std() ** 2
        )
    )
    return mean, covariance


def infer_parameters_from_mixture(expression_tensor, k=2):
    num_landmarks = expression_tensor.shape[0]
    device = expression_tensor.device

    results = Parallel(n_jobs=-1)(
        delayed(infer_landmark_parameters)(expression_tensor, i, k)
        for i in range(num_landmarks)
    )
    means, covariances = zip(*results)
    means = torch.stack(means).squeeze().to(device)
    covariances = torch.stack(covariances).squeeze().to(device)

    return means, covariances


def save_parameters(filepath, parameters_dict):
    with h5py.File(filepath, "w") as f:
        for landmark in parameters_dict.keys():
            group = f.create_group(landmark)

            write_mu(group, parameters_dict[landmark]["mu"])
            write_N(group, parameters_dict[landmark]["N"])

            S = group.create_group("S")

            write_S_off(S, parameters_dict[landmark]["S"]["off"])
            write_S_on(S, parameters_dict[landmark]["S"]["on"])


def load_parameters(filepath, landmarks):
    """Load parameters from an HDF5 file.

    Args:
        filepath (str): Path to the HDF5 file.

    Returns:
        tuple: Contains tensors for mu, N, and S.
    """
    with h5py.File(filepath, "r") as f:
        num_landmarks = len(landmarks)

        mu = torch.zeros((num_landmarks, 2))
        N = torch.zeros((num_landmarks, 2))
        S = torch.zeros((num_landmarks, 2, num_landmarks, num_landmarks))

        for i in range(num_landmarks):
            parameters = f[landmarks[i]]
            mu[i, :] = torch.tensor(parameters["mu"][()])
            N[i, :] = torch.tensor(parameters["N"][()])
            S[i, 0, :, :] = torch.tensor(parameters["S/off"][()])
            S[i, 1, :, :] = torch.tensor(parameters["S/on"][()])

    return mu, N, S


def calculate_bin_probability(
    expression_data: torch.Tensor,
    mean_values: torch.Tensor,
    covariances: torch.Tensor,
    landmark_indices: torch.Tensor,
    lambda_param: float,
) -> torch.Tensor:
    """
    Calculate the probability of a bin being in a particular state.

    Args:
        expression_data (torch.Tensor): Expression data (genes x samples).
        mean_values (torch.Tensor): Mean values for each landmark.
        covariances (torch.Tensor): Covariance values for each landmark.
        landmark_indices (torch.Tensor): Landmark indices for bin states.
        lambda_param (float): Regularization parameter for probabilities.

    Returns:
        torch.Tensor: Probabilities for the bin being in each state.
    """
    # Select mean and covariance values for the current bin
    mean_selected = mean_values.gather(1, landmark_indices.unsqueeze(1)).squeeze()
    covariance_selected = covariances.gather(1, landmark_indices.unsqueeze(1)).squeeze()

    # Construct a diagonal covariance matrix
    combined_covariance = torch.diag(covariance_selected)

    # Define multivariate normal distribution
    mvn = torch.distributions.MultivariateNormal(
        loc=mean_selected, covariance_matrix=combined_covariance
    )

    # Compute log probabilities for the bin
    log_probs = mvn.log_prob(expression_data.T)

    # Handle underflow or overflow in log probability scaling
    max_log_magnitude = torch.finfo(log_probs.dtype).max
    clipped_log_probs = torch.clamp(
        log_probs, min=-max_log_magnitude, max=max_log_magnitude
    )

    # Convert log probabilities to probabilities
    probabilities = torch.exp(clipped_log_probs)

    # Normalize probabilities to ensure they sum to 1
    if probabilities.sum() > 0:
        probabilities = probabilities / probabilities.sum()
        probabilities = kstat.utils.sparsegen(probabilities, lambda_param=lambda_param)
    else:
        probabilities = torch.zeros_like(probabilities)

    return probabilities


def calculate_bin_probabilities(
    expression_data: torch.Tensor,
    mean_values: torch.Tensor,
    covariances: torch.Tensor,
    unique_bins: torch.Tensor,
    lambda_param: float,
) -> torch.Tensor:
    """
    Calculate probabilities for each unique bin.

    Args:
        expression_data (torch.Tensor): Expression data (genes x samples).
        mean_values (torch.Tensor): Mean values for each landmark.
        covariances (torch.Tensor): Covariance values for each landmark.
        unique_bins (torch.Tensor): Unique bin indices (landmarks x bins).
        lambda_param (float): Regularization parameter for probabilities.

    Returns:
        torch.Tensor: Matrix of probabilities (samples x bins).
    """
    device = expression_data.device
    mean_values = mean_values.to(device)
    covariances = covariances.to(device)
    unique_bins = unique_bins.to(device)

    # Initialize probabilities matrix
    num_samples, num_bins = expression_data.shape[1], unique_bins.shape[1]
    bin_probabilities = torch.zeros((num_samples, num_bins), device=device)

    def process_bin(bin_idx):
        landmark_indices = unique_bins[:, bin_idx].long()
        probabilities = calculate_bin_probability(
            expression_data,
            mean_values,
            covariances,
            landmark_indices,
            lambda_param,
        )
        return bin_idx, probabilities

    # Use joblib to parallelize computation over bins
    results = Parallel(n_jobs=-1)(delayed(process_bin)(j) for j in range(num_bins))

    for bin_idx, probabilities in results:
        bin_probabilities[:, bin_idx] = probabilities

    # Set the first bin probabilities to zero (if required by the context)
    bin_probabilities[:, 0] = 0

    return bin_probabilities
