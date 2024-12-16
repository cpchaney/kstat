"""
Functions for modeling cell-cell communication.
"""

from typing import Tuple, List, Union
import torch
import pandas as pd
from scipy.sparse import coo_matrix
from geomloss import SamplesLoss
from skimage.filters import gaussian
import kstat.utils


def compute_mass(
    expression_tensor: torch.Tensor, genes: Union[str, List[str]], mask: coo_matrix
) -> torch.Tensor:
    """Compute mass of ligands or receptors.

    Args:
        expression_tensor (torch.Tensor): Tensor of expression data.
        genes (str or list): Gene(s) of interest.
        mask (coo_matrix): Mask for the expression tensor.

    Returns:
        torch.Tensor: Sparse tensor with applied mask and computed mass.
    """
    if isinstance(genes, str):
        filtered_expression = expression_tensor[genes].min(dim=0)[0]
    else:
        filtered_expression = expression_tensor[genes].min(dim=0)[0]
    masked = kstat.utils.apply_mask(filtered_expression, mask, return_dense=True)
    return torch.tensor(masked).to_sparse()


def calculate_interaction_loss(
    ligand_mass: torch.Tensor,
    receptor_mass: torch.Tensor,
    loss_fn: SamplesLoss,
    device: torch.device,
) -> float:
    """Calculate interaction loss between ligand and receptor masses.

    Args:
        ligand_mass (torch.Tensor): Ligand mass tensor.
        receptor_mass (torch.Tensor): Receptor mass tensor.
        loss_fn (SamplesLoss): Loss function from geomloss.
        device (torch.device): Device to move tensors to.

    Returns:
        float: Calculated interaction loss.
    """
    a = ligand_mass.values().to(device=device)
    x_i = ligand_mass.indices().T.float().to(device)
    b = receptor_mass.values().to(device)
    y_j = receptor_mass.indices().T.float().to(device)

    return loss_fn(a, x_i, b, y_j).item()


def get_ligand_receptor_interaction(
    interactor_table: pd.DataFrame,
    expression_data: Union[pd.DataFrame, torch.Tensor],
    X: torch.Tensor,
    X_r: torch.Tensor,
    mask: coo_matrix,
    device: torch.device,
    loss_fn: SamplesLoss,
) -> Tuple[List[float], List[float], List[float]]:
    """Calculate ligand-receptor interactions and their losses."""
    ligand_masses, receptor_masses, interaction_losses = [], [], []

    for interaction in interactor_table.index.unique():
        interaction_data = interactor_table.loc[interaction]

        receptor_genes = interaction_data[interaction_data["component"] == "receptor"][
            "gene"
        ]
        receptor_mass = compute_mass(X_r, receptor_genes, mask)
        receptor_masses.append(receptor_mass.sum().item())

        ligand_genes = interaction_data[interaction_data["component"] == "ligand"][
            "gene"
        ]
        ligand_mass = compute_mass(X, ligand_genes, mask)
        ligand_masses.append(ligand_mass.sum().item())

        interaction_loss = calculate_interaction_loss(
            ligand_mass, receptor_mass, loss_fn, device
        )
        interaction_losses.append(interaction_loss)

    return ligand_masses, receptor_masses, interaction_losses


def compute_receptive_field(
    cells: List[str],
    cell_type_table: pd.DataFrame,
    P: torch.Tensor,
    bin_keys: pd.DataFrame,
    height: int,
    width: int,
    tau: float,
) -> torch.Tensor:
    """Compute receptive field for given cells."""
    mask = torch.logical_not(
        torch.tensor(cell_type_table.index.isin(cells))[:, None].expand(-1, P.shape[1])
    )
    receptive_field = P.masked_fill(mask, 0).sum(dim=0)
    receptive_field = torch.where(receptive_field > 0, 1, 0)
    receptive_field = receptive_field[bin_keys.iloc[:, 1].to_numpy()]
    receptive_field = receptive_field.reshape(height, width)
    return torch.tensor(
        gaussian(receptive_field.float().numpy(), sigma=(tau / 4), preserve_range=True)
    )


def compute_ligand_masses(
    interactor_table: pd.DataFrame,
    expression_data: pd.DataFrame,
    X: torch.Tensor,
    mask: coo_matrix,
    receptive_field: torch.Tensor,
    height: int,
    width: int,
) -> List[float]:
    """Compute masses of ligands."""
    ligand_masses = []
    ligand_genes = interactor_table[interactor_table["component"] == "ligand"][
        "gene"
    ].unique()

    ligand_expectations = torch.zeros(len(ligand_genes), height * width)
    for i, ligand in enumerate(ligand_genes):
        ligand_expectations[i] = torch.tensor(
            kstat.utils.apply_mask(
                X[expression_data.index.get_loc(ligand)] * receptive_field,
                mask,
                return_dense=True,
            ).reshape(-1)
        )

    for interaction in interactor_table.index.unique():
        interaction_data = interactor_table.loc[interaction]
        ligand_genes = interaction_data[interaction_data["component"] == "ligand"][
            "gene"
        ]
        if isinstance(ligand_genes, str):
            mass = ligand_expectations[ligand_genes.index(ligand_genes)]
        else:
            indices = [ligand_genes.index(gene) for gene in ligand_genes]
            mass = ligand_expectations[indices].min(dim=0)[0]
        ligand_masses.append(mass.sum().item())

    return ligand_masses


def compute_interaction_scores(
    interactor_table: pd.DataFrame,
    ligand_masses: List[float],
    receptor_masses: List[float],
    interaction_losses: List[float],
) -> pd.DataFrame:
    """Compute interaction scores."""
    scores = (
        torch.tensor(ligand_masses) * torch.tensor(receptor_masses)
    ) / torch.tensor(interaction_losses)
    return pd.DataFrame(
        {
            "id_cp_interaction": interactor_table.index.unique(),
            "ligand_mass": ligand_masses,
            "receptor_mass": receptor_masses,
            "loss": interaction_losses,
            "score": scores.tolist(),
        }
    )
