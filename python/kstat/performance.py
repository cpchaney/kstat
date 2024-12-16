import geomloss
import numpy as np
import pandas as pd
from scipy.fft import fft2, ifft2
import scipy.io
import kstat.model
import kstat.utils
import torch
from typing import Tuple, List


def load_data(
    mask_path: str, landmark_measurements_path: str, landmarks_meta_path: str
) -> Tuple[torch.Tensor, torch.Tensor, pd.DataFrame]:
    """Load relevant data from the given paths."""
    mask = scipy.io.mmread(mask_path)
    landmark_measurements = scipy.io.mmread(landmark_measurements_path)
    L = kstat.utils.coo_to_torch(landmark_measurements, dense=False).coalesce()
    landmarks = pd.read_csv(landmarks_meta_path)
    return mask, L, landmarks


def configure_loss(p: int = 2, blur: float = 0.1) -> geomloss.SamplesLoss:
    """Configure the Sinkhorn divergence loss."""
    return geomloss.SamplesLoss("sinkhorn", p=p, blur=blur)


def compute_null_distributions(
    L: torch.Tensor,
    landmarks: pd.DataFrame,
    n_replicates: int,
    height: int,
    width: int,
    loss: geomloss.SamplesLoss,
    device: torch.device,
) -> pd.DataFrame:
    """Compute null distributions for landmark measurements."""
    N = torch.zeros(L.shape[0], n_replicates)

    for j in range(L.shape[0]):
        L_ = L[j].to_dense().reshape(height, width).to_sparse()
        y_j = L_.indices().T.float().to(device)

        random_samples = torch.randint(
            0, height, (n_replicates, y_j.shape[0], 2), device=device
        )
        for i, x_i in enumerate(random_samples):
            x_i = x_i.contiguous()
            N[j, i] = loss(x_i, y_j).cpu().item()

    return pd.DataFrame(data=N.numpy(), index=landmarks["gene"])


def calculate_reconstruction_error(
    X_i: torch.Tensor,
    L: torch.Tensor,
    loss: geomloss.SamplesLoss,
    unique_bins_key: List[int],
    height: int,
    width: int,
    device: torch.device,
) -> float:
    """Calculate reconstruction error between gene expression and landmarks."""
    X_i = X_i[unique_bins_key].reshape(height, width)
    X_i = kstat.utils.gaussian_filter(X_i)
    threshold = torch.sort(X_i.flatten(), descending=True)[0][int(L.sum().item()) - 1]
    X_i = torch.where(X_i >= threshold, X_i, 0)
    Y = L.to_dense().reshape(height, width).to_sparse()

    x_i = X_i.to_sparse().indices().T.float().to(device)
    y_j = Y.indices().T.float().to(device)

    return loss(x_i, y_j).cpu().item()


def calculate_reconstruction_errors(
    X: torch.Tensor,
    L: torch.Tensor,
    loss: geomloss.SamplesLoss,
    unique_bins_key: List[int],
    height: int,
    width: int,
    device: torch.device,
) -> torch.Tensor:
    """Calculate reconstruction errors for all cells."""
    return torch.tensor(
        [
            calculate_reconstruction_error(
                X[i], L, loss, unique_bins_key, height, width, device
            )
            for i in range(X.shape[0])
        ]
    )


def calculate_reconstruction_pvalues(
    null: pd.DataFrame, errors: torch.Tensor
) -> torch.Tensor:
    """Calculate p-values based on null distribution and reconstruction errors."""

    def compute_ecdf(data):
        sorted_data = np.sort(data)
        return lambda x: (sorted_data <= x).mean()

    ecdfs = [compute_ecdf(null.iloc[i].to_numpy()) for i in range(len(errors))]
    return torch.tensor([ecdf(error.item()) for ecdf, error in zip(ecdfs, errors)])


def calculate_loocv(
    expression_data: pd.DataFrame,
    common_landmarks: List[str],
    mu: torch.Tensor,
    N: torch.Tensor,
    S: torch.Tensor,
    L_: torch.sparse.Tensor,
    unique_bins: torch.Tensor,
    unique_bins_key: List[int],
    height: int,
    width: int,
    device: torch.device,
    loss: geomloss.SamplesLoss,
) -> torch.Tensor:
    """Perform Leave-One-Out Cross-Validation (LOOCV) for each landmark."""
    e = torch.zeros(L_.shape[0])

    for i, landmark in enumerate(common_landmarks):
        idx = np.delete(np.arange(len(common_landmarks)), i)
        E_ = torch.tensor(expression_data.loc[common_landmarks[idx]].to_numpy())
        mu_, N_, S_ = mu[idx], N[idx], S[idx, :, :][:, :, idx]
        P_ = kstat.model.compute_bin_probabilities(E_, mu_, N_, S_, unique_bins[idx])

        P_[:, 0] = 0  # Background probability adjustment
        X_i = torch.tensor(expression_data.loc[landmark].to_numpy()).unsqueeze(0) @ P_
        e[i] = calculate_reconstruction_error(
            X_i, L_, loss, unique_bins_key, height, width, device
        )

    return e


def calculate_spatial_correlation_fft(
    matrix1: np.ndarray, matrix2: np.ndarray
) -> np.ndarray:
    """Compute spatial correlation between two matrices using FFT."""
    if matrix1.shape != matrix2.shape:
        raise ValueError("The matrices must have the same dimensions")

    fft_matrix1 = fft2(matrix1)
    fft_matrix2 = fft2(matrix2)
    cross_power_spectrum = fft_matrix1 * np.conj(fft_matrix2)
    cross_correlation = ifft2(cross_power_spectrum).real

    norm = np.sqrt(np.sum(matrix1**2) * np.sum(matrix2**2))
    return cross_correlation / norm
