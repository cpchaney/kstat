# kstat/interactions/driver.py

from typing import Dict, List, Optional, Tuple

import pandas as pd
import torch

from kstat.interactions.expression import InteractionProjector
from kstat.interactions.transport import TransportCalculator


def make_labels(interaction_df: pd.DataFrame) -> Dict[int, str]:
    """
    Create human-readable labels for each interaction:
    e.g., "Lig1+Lig2-Rec1+Rec2".

    Args:
        interaction_df: DataFrame with component_type and mgi_symbol.

    Returns:
        Dictionary mapping id_cp_interaction -> label string.
    """
    labels = {}
    for iid, group in interaction_df.groupby("id_cp_interaction"):
        ligands = sorted(
            group.loc[group["component_type"] == "ligand", "mgi_symbol"].unique()
        )
        receptors = sorted(
            group.loc[group["component_type"] == "receptor", "mgi_symbol"].unique()
        )
        labels[iid] = "+".join(ligands) + "-" + "+".join(receptors)
    return labels


def run_interaction_pipeline(
    adata,
    interaction_df: pd.DataFrame,
    cell_type_list: List[str],
    cell_autonomous: bool,
    cells_bins_prob: torch.Tensor,
    unique_bins_key: torch.Tensor,
    height: int,
    width: int,
    device: str = "cpu",
    p: int = 1,
    blur: float = 0.05,
    reach: Optional[float] = None,
    scaling: float = 0.9,
    mode: str = "balanced",  # or "unbalanced"
    normalize_before_unbalanced: Optional[bool] = None,
    k: float = 1.0,
    return_plans: bool = False,
) -> Tuple[pd.DataFrame, Dict[int, torch.Tensor]]:
    """
    Full pipeline to score spatial ligand-receptor interactions using optimal transport.

    Projects expression to spatial bins and computes transport-based similarity between
    ligand and receptor distributions.

    Args:
        adata: AnnData object with expression data.
        interaction_df: DataFrame of interaction annotations.
        cell_type_list: Receptor cell types.
        cell_autonomous: If False, ligand cells exclude receptor cells.
        cells_bins_prob: (n_cells, n_bins) cell-to-bin probability matrix.
        unique_bins_key: Maps spatial bins to compact index set (long tensor).
        height, width: Image/grid dimensions.
        device: Torch device.
        p: OT cost exponent (e.g., 1 or 2).
        blur: Entropic regularization strength.
        reach: Optional reach parameter for unbalanced OT.
        scaling: GeomLoss scaling factor (controls convergence).
        mode: "balanced" or "unbalanced".
        normalize_before_unbalanced: Optional manual override.
        k: Score scaling factor.
        return_plans: Whether to return full OT transport maps.

    Returns:
        results_df: DataFrame with interaction metrics.
        plans_dict: Optional, contains transport plans and raw mass maps.
    """
    device_t = torch.device(device)

    # Set normalization logic for unbalanced mode
    if normalize_before_unbalanced is None:
        normalize_before_unbalanced = mode == "unbalanced"

    # 1. Setup projector + OT solver
    projector = InteractionProjector(adata, interaction_df)
    calculator = TransportCalculator(
        unique_bins_key=unique_bins_key,
        height=height,
        width=width,
        device=device_t,
        p=p,
        blur=blur,
        reach=(reach if mode == "unbalanced" else None),
        scaling=scaling,
        mode=mode,
        normalize_before_unbalanced=normalize_before_unbalanced,
    )

    # 2. Project ligand/receptor expression to spatial bins
    projected = projector.project_all_interactions(
        cells_bins_prob=cells_bins_prob,
        receiving_cell_types=cell_type_list,
        cell_autonomous=cell_autonomous,
    )

    results: List[Dict] = []
    plans_dict: Dict = {}
    eps = 1e-12  # Small constant to prevent division by zero

    # 3. Score each interaction using OT
    for interaction_id, (lig_expr, rec_expr) in projected.items():
        lig_mass_raw = float(lig_expr.sum().item())
        rec_mass_raw = float(rec_expr.sum().item())

        if lig_mass_raw == 0.0 or rec_mass_raw == 0.0:
            continue  # Skip interactions with no expression

        (
            loss,
            ligand_full,
            receptor_full,
            delivered_prop,
            received_prop,
        ) = calculator.compute_loss_and_plan(
            ligand_expr=lig_expr,
            receptor_expr=rec_expr,
        )

        if not (loss == loss):  # NaN guard
            continue

        # Coupling-aware interaction score (inversely proportional to cost)
        score = k * torch.sqrt(ligand_full.sum() * receptor_full.sum()) / (loss + eps)

        results.append(
            {
                "id_cp_interaction": interaction_id,
                "ligand_mass_full": ligand_full.sum().item(),
                "receptor_mass_full": receptor_full.sum().item(),
                "loss": loss,
                "score": score.item(),
            }
        )

        if return_plans:
            plans_dict[interaction_id] = {
                "ligand": ligand_full.cpu(),
                "receptor": receptor_full.cpu(),
                "delivered_prop": delivered_prop.cpu(),
                "received_prop": received_prop.cpu(),
            }

    # 4. Build DataFrame with results and add labels
    results_df = pd.DataFrame(results)

    labels_map = make_labels(interaction_df)
    if not results_df.empty:
        results_df["label"] = results_df["id_cp_interaction"].map(labels_map)
    else:
        results_df["label"] = []

    return results_df, plans_dict
