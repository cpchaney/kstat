# tests/test_expression.py

import anndata as ad
import numpy as np
import pandas as pd
import torch

from kstat.interactions.expression import (
    InteractionProjector,
    mask_cells_by_type,
    project_expression_to_bins,
)


def test_mask_cells_by_type():
    """
    Test that mask_cells_by_type correctly identifies matching cell types.
    """
    # Create dummy AnnData object with 4 cells and a cell type annotation
    adata = ad.AnnData(np.zeros((4, 2)))
    adata.obs["cell_type_short"] = ["A", "B", "A", "C"]

    # Apply mask for cell type "A"
    mask = mask_cells_by_type(adata, ["A"])

    # Expected: Only rows 0 and 2 are type A
    assert mask.tolist() == [True, False, True, False]


def test_project_expression_to_bins():
    """
    Test projection of expression data to spatial bins using a binary cell mask.
    """
    # Dummy gene expression matrix (3 cells × 2 genes)
    expr = np.array([[1, 2], [3, 4], [5, 6]])

    # Mask: select cells 0 and 2
    mask = np.array([True, False, True])

    # Dummy cell-to-bin probabilities (equal weight to all 4 bins)
    prob = np.ones((3, 4)) / 3

    # Test projecting both genes to bins
    gene_indices = [0, 1]
    out = project_expression_to_bins(expr, mask, prob, gene_indices, reduce="min")

    assert isinstance(out, torch.Tensor)
    assert out.shape == (4,)
