from typing import Literal, Optional, Tuple, Union

import torch
from geomloss import SamplesLoss
from pykeops.torch import LazyTensor

# Type hint for the OT mode
Mode = Literal["balanced", "unbalanced"]


class TransportCalculator:
    """
    Computes optimal transport (OT)-based interactions (e.g. ligand-receptor communication)
    between two spatial signals using Sinkhorn divergence and dual potentials via KeOps.
    """

    def __init__(
        self,
        unique_bins_key: torch.Tensor,
        height: int,
        width: int,
        device: Union[str, torch.device] = "cpu",
        p: int = 1,
        blur: float = 0.05,
        reach: Optional[float] = None,
        scaling: float = 0.9,
        mode: Mode = "balanced",
        normalize_before_unbalanced: bool = True,
        eps_zero: float = 0.0,
        dtype: torch.dtype = torch.float32,
    ):
        """
        Initialize the transport solver with OT parameters and grid layout.

        Args:
            unique_bins_key (Tensor): Maps full grid (H*W) to compact vector indices, or -1.
            height, width (int): Spatial grid dimensions.
            device (str or torch.device): Device for computation.
            p (int): OT power (1 or 2).
            blur (float): Blur parameter for Sinkhorn kernel.
            reach (float or None): Unbalanced transport penalty.
            scaling (float): Scaling factor for Sinkhorn iteration.
            mode (str): 'balanced' or 'unbalanced' OT.
            normalize_before_unbalanced (bool): Normalize before unbalanced transport.
            eps_zero (float): Drop values below this threshold.
            dtype (torch.dtype): Floating point precision.
        """
        self.device = torch.device(device)
        self.dtype = dtype
        self.height, self.width = height, width
        self.n_total_bins = height * width

        self.unique_bins_key = unique_bins_key.to(self.device).long()

        self.p = int(p)
        self.blur = float(blur)
        self.scaling = float(scaling)
        self.mode = mode
        self.normalize_before_unbalanced = bool(normalize_before_unbalanced)
        self.eps_zero = float(eps_zero)

        # Validate mode and reach
        if self.mode == "unbalanced":
            if reach is None or not torch.isfinite(torch.tensor(reach)):
                raise ValueError("Provide a finite 'reach' for unbalanced OT.")
        self.reach_for_solver = None if self.mode == "balanced" else reach

        # Coordinates of all grid positions (flattened)
        self.all_coords = self._make_coords(height, width).to(self.device, dtype=self.dtype)

        # OT loss functions (with and without duals)
        self.loss_fn = SamplesLoss(
            loss="sinkhorn",
            p=self.p,
            blur=self.blur,
            reach=self.reach_for_solver,
            scaling=self.scaling,
            debias=False,
        )
        self._pot_solver = SamplesLoss(
            loss="sinkhorn",
            p=self.p,
            blur=self.blur,
            reach=self.reach_for_solver,
            scaling=self.scaling,
            debias=False,
            potentials=True,
        )

    # ---------- Helper Methods ----------

    def _make_coords(self, h: int, w: int) -> torch.Tensor:
        """Generate (x, y) coordinates for the HxW grid as a (H*W, 2) tensor."""
        y, x = torch.meshgrid(
            torch.arange(h, dtype=self.dtype),
            torch.arange(w, dtype=self.dtype),
            indexing="ij"
        )
        return torch.stack([x.reshape(-1), y.reshape(-1)], dim=-1)

    def _expand_to_full_space(self, compact: torch.Tensor) -> torch.Tensor:
        """
        Expand a compact vector (K,) to the full grid (H*W,) using unique_bins_key.
        """
        compact = compact.to(self.device, dtype=self.dtype)
        full = torch.zeros(self.n_total_bins, device=self.device, dtype=self.dtype)
        valid = self.unique_bins_key >= 0
        full[valid] = compact[self.unique_bins_key[valid]]
        return full

    def _sparsify(self, weights: torch.Tensor):
        """
        Extract nonzero values and their coordinates for OT computation.
        Returns (weights, coords, indices).
        """
        mask = (weights > self.eps_zero) if self.eps_zero > 0 else (weights > 0)
        idx = torch.nonzero(mask, as_tuple=False).flatten()
        if idx.numel() == 0:
            return weights.new_zeros((0,)), self.all_coords.new_zeros((0, 2)), idx
        return weights[idx], self.all_coords[idx], idx

    def _compute_transported_prop(
        self,
        coords: torch.Tensor,
        transported: torch.Tensor,
        signaller_full: torch.Tensor
    ) -> torch.Tensor:
        """
        Map transported mass back to grid and normalize it relative to signaller distribution.
        Avoids NaNs by safe division and clamps output to [0,1].
        """
        x_int = coords[:, 0].long()
        y_int = coords[:, 1].long()
        flat_indices = (y_int * self.width + x_int).to(self.device)

        flat = torch.zeros(self.n_total_bins, dtype=transported.dtype, device=self.device)
        flat.scatter_add_(0, flat_indices, transported)

        sig = signaller_full.to(device=self.device, dtype=transported.dtype)
        sig_sum = sig.sum()
        if sig_sum.item() == 0.0:
            return torch.zeros_like(flat)

        sig = sig / sig_sum
        mask = sig > 0
        out = torch.zeros_like(flat)
        out[mask] = (flat[mask] / sig[mask]).clamp_(max=1.0)
        return out

    # ---------- Main Method ----------

    @torch.no_grad()
    def compute_loss_and_plan(
        self,
        ligand_expr: torch.Tensor,
        receptor_expr: torch.Tensor,
        normalize_balanced: bool = True
    ) -> Tuple[float, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Compute the OT loss and transport plan between two spatial signals.

        Args:
            ligand_expr (Tensor): Compact vector of ligand signal (K,).
            receptor_expr (Tensor): Compact vector of receptor signal (K,).
            normalize_balanced (bool): Ignored, included for API compatibility.

        Returns:
            Tuple containing:
                - loss (float)
                - ligand_full (H*W,)
                - receptor_full (H*W,)
                - delivered_prop (H*W,)
                - received_prop (H*W,)
        """
        # Expand to full grid
        ligand_full = self._expand_to_full_space(ligand_expr)
        receptor_full = self._expand_to_full_space(receptor_expr)

        # Sparsify inputs
        a_raw, x_coords, idx_a = self._sparsify(ligand_full)
        b_raw, y_coords, idx_b = self._sparsify(receptor_full)

        # Early return on empty support
        if (
            a_raw.numel() == 0 or b_raw.numel() == 0
            or a_raw.sum().item() == 0.0
            or b_raw.sum().item() == 0.0
        ):
            empty = ligand_full.new_zeros((self.n_total_bins,))
            return float("nan"), ligand_full, receptor_full, empty, empty

        # Normalize (balanced or unbalanced)
        if self.mode == "balanced" or self.normalize_before_unbalanced:
            a = a_raw / a_raw.sum()
            b = b_raw / b_raw.sum()
        else:
            a, b = a_raw, b_raw

        # Prepare inputs for KeOps
        x = x_coords.to(self.device, dtype=self.dtype).contiguous()
        y = y_coords.to(self.device, dtype=self.dtype).contiguous()
        a = a.to(self.device, dtype=self.dtype).contiguous()
        b = b.to(self.device, dtype=self.dtype).contiguous()

        # Compute Sinkhorn divergence (OT loss)
        loss = self.loss_fn(a, x, b, y).item()

        # Get dual potentials
        F, G = self._pot_solver(a, x, b, y)
        F, G = F.to(self.device), G.to(self.device)

        # Define distance matrix using LazyTensor
        x_i = LazyTensor(x[:, None, :])
        y_j = LazyTensor(y[None, :, :])
        sq = ((x_i - y_j) ** 2).sum(-1)

        # Cost function
        if self.p == 1:
            C_ij = sq.sqrt()
        elif self.p == 2:
            C_ij = 0.5 * sq
        else:
            C_ij = (sq ** (self.p / 2.0)) / self.p

        # Compute transport plan
        eps = self.blur ** self.p
        plan = ((F.view(-1, 1) + G.view(1, -1) - C_ij) / eps).exp()
        plan *= a.view(-1, 1) * b.view(1, -1)

        delivered = plan.sum(dim=1)
        received = plan.sum(dim=0)

        # Map transported mass back to grid
        delivered_prop = self._compute_transported_prop(x, delivered, ligand_full)
        received_prop = self._compute_transported_prop(y, received, receptor_full)

        return loss, ligand_full, receptor_full, delivered_prop, received_prop
