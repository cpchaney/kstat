# tests/test_transport.py

import torch
from kstat.interactions.transport import TransportCalculator


def test_transport_loss():
    """
    Unit test for the compute_loss_and_plan method in TransportCalculator.

    Sets up a minimal toy grid with ligand and receptor expressions,
    then checks that the OT loss is computed and is a non-negative float.
    """
    # Define mapping from full grid to compact representation (-1 would mean "no data")
    unique_bins = torch.tensor([0, 1, 2, 3])

    # Define spatial grid size (2 rows × 3 cols)
    height, width = 2, 3

    # Define compact ligand and receptor signals (length = 4)
    lig = torch.tensor([0.1, 0.3, 0.0, 0.6], dtype=torch.float32)
    rec = torch.tensor([0.0, 0.2, 0.4, 0.4], dtype=torch.float32)

    # Initialize the TransportCalculator with default Sinkhorn OT settings
    tc = TransportCalculator(unique_bins, height, width)

    # Compute OT-based loss and other outputs
    result = tc.compute_loss_and_plan(lig, rec)

    # Ensure the output is a 5-tuple: (loss, ligand_full, receptor_full, delivered, received)
    assert isinstance(result, tuple)
    assert len(result) == 5

    loss = result[0]

    # Validate the returned loss
    assert isinstance(loss, float)
    assert loss >= 0.0
