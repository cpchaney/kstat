from typing import Literal, Optional, Tuple, Union

import torch
from geomloss import SamplesLoss
from pykeops.torch import LazyTensor

Mode = Literal["balanced", "unbalanced"]


class TransportCalculator:
    def __init__(
        self,
        unique_bins_key: torch.Tensor,  # (H*W,) int, -1 means "no data", else index into compact vectors
        height: int,
        width: int,
        device: Union[str, torch.device] = "cpu",
        p: int = 1,
        blur: float = 0.05,
        reach: Optional[float] = None,  # only used in unbalanced mode
        scaling: float = 0.9,
        mode: Mode = "balanced",
        normalize_before_unbalanced: bool = True,
        eps_zero: float = 0.0,  # >0 to drop tiny weights
        dtype: torch.dtype = torch.float32,
    ):
        self.device = torch.device(device)
        self.dtype = dtype
        self.height, self.width = height, width
        self.n_total_bins = height * width

        # unique_bins_key maps full grid -> compact vector idx (or -1)
        self.unique_bins_key = unique_bins_key.to(self.device).long()

        self.p = int(p)
        self.blur = float(blur)
        self.scaling = float(scaling)
        if mode not in ("balanced", "unbalanced"):
            raise ValueError("mode must be 'balanced' or 'unbalanced'")
        self.mode = mode
        self.normalize_before_unbalanced = bool(normalize_before_unbalanced)
        self.eps_zero = float(eps_zero)

        # Full-grid coordinates (float for OT)
        self.all_coords = self._make_coords(height, width).to(
            self.device, dtype=self.dtype
        )

        # Unbalanced config
        self.reach_for_solver = None if mode == "balanced" else reach
        if mode == "unbalanced":
            if reach is None or not torch.isfinite(torch.tensor(reach)):
                raise ValueError("Provide a finite 'reach' for unbalanced OT.")

        # GeomLoss solvers
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

    # ---------- helpers ----------
    def _make_coords(self, h: int, w: int) -> torch.Tensor:
        y, x = torch.meshgrid(
            torch.arange(h, dtype=self.dtype),
            torch.arange(w, dtype=self.dtype),
            indexing="ij",
        )
        return torch.stack([x.reshape(-1), y.reshape(-1)], dim=-1)  # (H*W, 2)

    def _expand_to_full_space(self, compact: torch.Tensor) -> torch.Tensor:
        """Map a length-K compact vector to a length H*W flat vector via unique_bins_key."""
        compact = compact.to(self.device, dtype=self.dtype)
        full = torch.zeros(self.n_total_bins, device=self.device, dtype=self.dtype)
        valid = self.unique_bins_key >= 0
        full[valid] = compact[self.unique_bins_key[valid]]
        return full

    def _sparsify(self, weights: torch.Tensor):
        """Return (w, coords, idx) where w>eps and idx are flat indices in the full grid."""
        w = weights
        mask = (w > self.eps_zero) if self.eps_zero > 0.0 else (w > 0)
        idx = torch.nonzero(mask, as_tuple=False).flatten()
        if idx.numel() == 0:
            # empty
            return w.new_zeros((0,)), self.all_coords.new_zeros((0, 2)), idx
        return w[idx], self.all_coords[idx], idx

    def _compute_transported_prop(
        self,
        coords: torch.Tensor,  # (N,2) float coords exactly on integer grid
        transported: torch.Tensor,  # (N,)
        signaller_full: torch.Tensor,  # (H*W,)
    ) -> torch.Tensor:
        """Scatter-add transported mass to full grid, divide by normalized signaller, clamp [0,1]."""
        device = transported.device
        dtype = transported.dtype

        # coords are integer-valued stored as float; cast safely
        x_int = coords[:, 0].long()
        y_int = coords[:, 1].long()
        flat_indices = (y_int * self.width + x_int).to(device)

        # scatter-add to allow multiple hits per bin
        flat = torch.zeros(self.n_total_bins, dtype=dtype, device=device)
        flat.scatter_add_(0, flat_indices, transported)

        sig = signaller_full.to(device=device, dtype=dtype)
        sig_sum = sig.sum()
        if sig_sum.item() == 0.0:
            # No signal anywhere -> return zeros, avoids NaNs
            return torch.zeros_like(flat)

        sig = sig / sig_sum
        mask = sig > 0
        out = torch.zeros_like(flat)
        out[mask] = (flat[mask] / sig[mask]).clamp_(max=1.0)
        return out

    # ---------- main ----------
    @torch.no_grad()
    def compute_loss_and_plan(
        self,
        ligand_expr: torch.Tensor,  # (K_ligand,)
        receptor_expr: torch.Tensor,  # (K_receptor,)
        normalize_balanced: bool = True,  # kept for API compat
    ) -> Tuple[float, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Returns:
            loss (float),
            ligand_full (H*W,),
            receptor_full (H*W,),
            delivered_prop (H*W,),
            received_prop (H*W,)
        """
        # 1) Expand compact vectors to full grid
        ligand_full = self._expand_to_full_space(ligand_expr)
        receptor_full = self._expand_to_full_space(receptor_expr)

        # 2) Sparsify
        a_raw, x_coords, idx_a = self._sparsify(ligand_full)
        b_raw, y_coords, idx_b = self._sparsify(receptor_full)

        # 3) Handle empty supports consistently
        if (
            (a_raw.numel() == 0)
            or (b_raw.numel() == 0)
            or (a_raw.sum().item() == 0.0)
            or (b_raw.sum().item() == 0.0)
        ):
            empty = ligand_full.new_zeros((self.n_total_bins,))
            return float("nan"), ligand_full, receptor_full, empty, empty

        # 4) Normalize masses according to mode
        if self.mode == "balanced" or self.normalize_before_unbalanced:
            a = a_raw / a_raw.sum()
            b = b_raw / b_raw.sum()
        else:
            a, b = a_raw, b_raw

        # 5) Prepare tensors for KeOps
        x = x_coords.to(self.device, dtype=self.dtype).contiguous()
        y = y_coords.to(self.device, dtype=self.dtype).contiguous()
        a = a.to(self.device, dtype=self.dtype).contiguous()
        b = b.to(self.device, dtype=self.dtype).contiguous()

        # 6) OT loss
        loss = self.loss_fn(a, x, b, y).item()

        # 7) Dual potentials
        F, G = self._pot_solver(a, x, b, y)
        F = F.to(self.device, dtype=self.dtype).contiguous()
        G = G.to(self.device, dtype=self.dtype).contiguous()

        N, M = x.shape[0], y.shape[0]
        eps = self.blur**self.p

        # 8) Kernelized plan proxy via duals
        x_i = LazyTensor(x[:, None, :])  # (N,1,D)
        y_j = LazyTensor(y[None, :, :])  # (1,M,D)
        sq = ((x_i - y_j) ** 2).sum(-1)

        if self.p == 1:
            C_ij = sq.sqrt()
        elif self.p == 2:
            C_ij = 0.5 * sq
        else:
            C_ij = (sq ** (self.p / 2.0)) / self.p

        F_i = LazyTensor(F.view(N, 1, 1))
        G_j = LazyTensor(G.view(1, M, 1))
        a_i = LazyTensor(a.view(N, 1, 1))
        b_j = LazyTensor(b.view(1, M, 1))

        plan = ((F_i + G_j - C_ij) / eps).exp() * (a_i * b_j)  # (N,M,1)

        delivered = plan.sum(dim=1).squeeze(-1)  # (N,)
        received = plan.sum(dim=0).squeeze(-1)  # (M,)

        # 9) Proportions on full grid (scatter-add + safe ratio)
        delivered_prop = self._compute_transported_prop(x, delivered, ligand_full)
        received_prop = self._compute_transported_prop(y, received, receptor_full)

        return (
            loss,
            ligand_full,
            receptor_full,
            delivered_prop,
            received_prop,
        )
