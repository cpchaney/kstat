from typing import Dict, List, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy.sparse import issparse
from tqdm import tqdm

# ----------------------------- #
#         HELPER FUNCTIONS     #
# ----------------------------- #

def find_signal_knee(signal: torch.Tensor) -> float:
    """
    Detects the 'knee' (elbow point) of a 1D signal vector using
    the maximum perpendicular distance method.

    Args:
        signal (torch.Tensor): 1D vector of signal values.

    Returns:
        float: Threshold value at the detected knee point.
    """
    signal_unique_sorted = torch.unique(signal).sort()[0]

    x = torch.arange(len(signal_unique_sorted), dtype=signal.dtype, device=signal.device)
    y = signal_unique_sorted

    point1 = torch.stack([x[0], y[0]])
    point2 = torch.stack([x[-1], y[-1]])
    line_vec = point2 - point1
    line_vec_norm = line_vec / torch.norm(line_vec)

    vecs_to_line = torch.stack([x - point1[0], y - point1[1]], dim=1)
    proj_lengths = (vecs_to_line * line_vec_norm).sum(dim=1)
    proj_points = point1.unsqueeze(0) + proj_lengths.unsqueeze(1) * line_vec_norm.unsqueeze(0)

    actual_points = torch.stack([x, y], dim=1)
    distances = torch.norm(actual_points - proj_points, dim=1)

    knee_index = torch.argmax(distances).item()
    return y[knee_index].item()


def mask_cells_by_type(
    adata: ad.AnnData, cell_types: Union[str, List[str]]
) -> np.ndarray:
    """
    Return a boolean mask for selecting cells by type.

    Args:
        adata (AnnData): Annotated single-cell dataset.
        cell_types (str or list): Cell type(s) to select.

    Returns:
        np.ndarray: Boolean mask of shape (n_cells,).
    """
    if isinstance(cell_types, str):
        cell_types = [cell_types]
    return adata.obs["cell_type_short"].isin(cell_types).to_numpy()


def project_expression_to_bins(
    expr_matrix: Union[np.ndarray, torch.Tensor],
    cell_mask: np.ndarray,
    cell_bin_probs: Union[np.ndarray, torch.Tensor],
    gene_indices: Union[List[int], np.ndarray],
    reduce: str = "min",
) -> torch.Tensor:
    """
    Project gene expression from selected cells to spatial bins.

    Args:
        expr_matrix: (n_cells, n_genes) gene expression values.
        cell_mask: (n_cells,) boolean mask selecting cell types.
        cell_bin_probs: (n_cells, n_bins) cell-to-bin probabilities.
        gene_indices: Indices of selected genes.
        reduce: 'min' (default), 'sum', or 'mean' to combine genes.

    Returns:
        (n_bins,) torch.Tensor of projected expression.
    """
    if isinstance(expr_matrix, np.ndarray):
        expr_matrix = torch.tensor(expr_matrix, dtype=torch.float32)
    if isinstance(cell_bin_probs, np.ndarray):
        cell_bin_probs = torch.tensor(cell_bin_probs, dtype=torch.float32)

    expr_subset = expr_matrix[cell_mask][:, gene_indices]
    probs_subset = cell_bin_probs[cell_mask]

    expected_expr = torch.matmul(expr_subset.T, probs_subset)  # (n_genes, n_bins)

    if reduce == "min":
        return expected_expr.min(dim=0).values
    elif reduce == "sum":
        return expected_expr.sum(dim=0)
    elif reduce == "mean":
        return expected_expr.mean(dim=0)
    else:
        raise ValueError(f"Unsupported reduce mode: {reduce}")


# ----------------------------- #
#      MAIN PROJECTOR CLASS    #
# ----------------------------- #

class InteractionProjector:
    """
    Projects ligand-receptor gene expression into spatial bins
    using probabilistic cell location estimates.
    """

    def __init__(
        self,
        adata: ad.AnnData,
        interaction_df: pd.DataFrame,
        gene_key: str = "mgi_symbol",
    ):
        """
        Args:
            adata (AnnData): Single-cell expression dataset.
            interaction_df (pd.DataFrame): Table with columns:
                ['id_cp_interaction', 'mgi_symbol', 'component_type'].
            gene_key (str): Key to match genes (default: 'mgi_symbol').
        """
        self.adata = adata
        self.interaction_df = interaction_df
        self.gene_index = {g: i for i, g in enumerate(adata.var_names)}

    def _get_gene_indices(self, gene_names: List[str]) -> List[int]:
        """Convert gene names to indices in adata.var_names."""
        return [self.gene_index[g] for g in gene_names if g in self.gene_index]

    def project_all_interactions(
        self,
        cells_bins_prob: Union[np.ndarray, torch.Tensor],
        receiving_cell_types: Union[str, List[str]],
        cell_autonomous: bool = True,
    ) -> Dict[int, Tuple[torch.Tensor, torch.Tensor]]:
        """
        Project expression of all ligand-receptor interactions to spatial bins.

        Args:
            cells_bins_prob: (n_cells, n_bins) cell-to-bin probability matrix.
            receiving_cell_types (str or list): Cells expressing receptors.
            cell_autonomous (bool): If True, ligands can come from all cells;
                                    if False, only from non-receptor cells.

        Returns:
            Dict[int, Tuple[ligand_tensor, receptor_tensor]]
        """
        # Load gene expression matrix
        expr_matrix = self.adata.X.toarray() if issparse(self.adata.X) else self.adata.X

        # Define masks
        receptor_mask = mask_cells_by_type(self.adata, receiving_cell_types)
        ligand_mask = np.ones(expr_matrix.shape[0], dtype=bool) if cell_autonomous else ~receptor_mask

        result = {}
        grouped = self.interaction_df.groupby("id_cp_interaction")

        for interaction_id, group in tqdm(grouped, desc="Projecting expression"):
            ligands = group.loc[group["component_type"] == "ligand", "mgi_symbol"]
            receptors = group.loc[group["component_type"] == "receptor", "mgi_symbol"]

            ligand_idx = self._get_gene_indices(ligands)
            receptor_idx = self._get_gene_indices(receptors)

            if not ligand_idx or not receptor_idx:
                continue  # Skip interactions with missing genes

            # Project ligand expression
            ligand_expr = project_expression_to_bins(
                expr_matrix, ligand_mask, cells_bins_prob, ligand_idx, reduce="min"
            )
            ligand_threshold = find_signal_knee(ligand_expr)
            ligand_expr[ligand_expr < ligand_threshold] = 0

            # Project receptor expression
            receptor_expr = project_expression_to_bins(
                expr_matrix, receptor_mask, cells_bins_prob, receptor_idx, reduce="min"
            )
            receptor_threshold = find_signal_knee(receptor_expr)
            receptor_expr[receptor_expr < receptor_threshold] = 0

            result[interaction_id] = (ligand_expr, receptor_expr)

        return result
