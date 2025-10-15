# mouse_e18.5d_kidney_ligand_receptor_inference.py

"""
Run spatial ligand-receptor interaction inference on mouse E18.5d kidney data.
Performs interaction projection, Sinkhorn transport, and CellPhoneDB compatibility export.

Usage:
    python mouse_e18.5d_kidney_ligand_receptor_inference.py
"""

import os
import pickle
import sqlite3
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import yaml
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
from scipy.io import mmread

from kstat.interactions.driver import run_interaction_pipeline

# -------------------- #
# Save utilities       #
# -------------------- #

def save_distributions_to_hdf5(data_dict, file_path):
    """Save transport distributions (as tensors) to HDF5 format."""
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with h5py.File(file_path, "w") as f:
        for interaction_id, interaction_data in data_dict.items():
            grp = f.create_group(interaction_id)
            for key, tensor in interaction_data.items():
                if not isinstance(tensor, torch.Tensor):
                    raise TypeError(f"{interaction_id}/{key} is not a torch.Tensor")
                grp.create_dataset(key, data=tensor.detach().cpu().numpy())

# -------------------- #
# Load configuration   #
# -------------------- #

with open("../../config/config.yaml") as f:
    config = yaml.safe_load(f)

PROJECT_ROOT = Path(config["project_root"])
NGS_ROOT = Path(config["ngs_root"])

# -------------------- #
# Load interaction DB  #
# -------------------- #

sql_path = PROJECT_ROOT / config["sql_query_path"]
with open(sql_path, "r") as f:
    sql_query = f.read()

db_path = Path.home() / "databases" / config["database"]
with sqlite3.connect(db_path) as conn:
    interaction_df = pd.read_sql_query(sql_query, conn)

print(f"Loaded {len(interaction_df)} ligand-receptor pairs.")
interaction_df.to_csv(PROJECT_ROOT / "data/ligand_receptor_components.csv", index=False)

# ------------------------------ #
# Load cell-bin probabilities    #
# ------------------------------ #

cells_bins_probs = mmread(PROJECT_ROOT / "data/mouse_e18.5d_cells_bins_probabilities.mtx").tocsr()
cells_bins_probs = torch.tensor(cells_bins_probs.toarray(), dtype=torch.float32)

# ------------------------------ #
# Load AnnData and filter genes  #
# ------------------------------ #

adata = sc.read_h5ad(PROJECT_ROOT / "../kidney_development/mouse_e18.5d/data/mouse_e18.5d_sn.h5ad")
cell_type_table = pd.read_csv(
    NGS_ROOT / "archive/STAT/data/foxd1-cre_myc_mycn_wild_type_cell_types.csv",
    index_col=0,
)
adata.obs["cell_type_short"] = pd.Categorical(
    cell_type_table.loc[adata.obs_names]["cell.type.short"]
)

# Keep only genes in interaction list
mgi_genes = pd.unique(interaction_df["mgi_symbol"].values)
adata = adata[:, adata.var_names.isin(mgi_genes)].copy()

# Set appropriate layers
adata.layers["counts"] = adata.X.copy()
adata.X = adata.layers["imputed"]  # Replace this if another expression layer is preferred

# -------------------- #
# Load bin mapping     #
# -------------------- #

unique_bins_key_table = pd.read_csv(PROJECT_ROOT / "data/mouse_e18.5d_unique_bins_key.csv")
unique_bins_key = torch.tensor(unique_bins_key_table["index"])

# -------------------- #
# Run interaction pipeline
# -------------------- #

cell_type_list = ["NPC1", "NPC2", "NPC3"]  # Receiving cell types
device = "cuda" if torch.cuda.is_available() else "cpu"

start_time = time.time()
results_df, plans = run_interaction_pipeline(
    adata=adata,
    interaction_df=interaction_df,
    cell_type_list=cell_type_list,
    cell_autonomous=False,
    cells_bins_prob=cells_bins_probs,
    unique_bins_key=unique_bins_key,
    height=config["scaled_height"],
    width=config["scaled_width"],
    device=device,
    p=config.get("p", 1),
    blur=config.get("blur", 3.0),
    reach=config.get("reach", 20.0),
    scaling=config.get("scaling"),
    mode="unbalanced",
    normalize_before_unbalanced=True,
    return_plans=True,
)
elapsed_time = time.time() - start_time
print(f"Pipeline runtime: {elapsed_time:.2f} seconds")

# -------------------- #
# Save pipeline results
# -------------------- #

cell_type_str = "_".join(cell_type.lower() for cell_type in cell_type_list)

results_path = PROJECT_ROOT / f"output/mouse_e18.5d_kidney_{cell_type_str}_interactions.csv"
results_df.to_csv(results_path, index=False)
print(f"Saved results to {results_path}")

# Save transport distributions to both pickle and HDF5
pkl_path = PROJECT_ROOT / f"data/mouse_e18.5d_kidney_{cell_type_str}_receptor_saturations.pkl"
with open(pkl_path, "wb") as f:
    pickle.dump(plans, f, protocol=pickle.HIGHEST_PROTOCOL)

h5_path = PROJECT_ROOT / f"data/mouse_e18.5d_kidney_{cell_type_str}_transport_distributions.h5"
save_distributions_to_hdf5(plans, str(h5_path))

# -------------------- #
# Prepare CellPhoneDB  #
# -------------------- #

# Tag receiver/sender roles
adata.obs["role"] = "sender"
adata.obs.loc[adata.obs["cell_type_short"].isin(cell_type_list), "role"] = "receiver"
adata.X = adata.layers["counts"].copy()
adata.raw = sc.AnnData(
    X=adata.layers["counts"].copy(),
    var=adata.var.copy(),
    obs=adata.obs.copy(),
)

# Export metadata
metadata_df = pd.DataFrame({
    "Cell": adata.obs_names,
    "cell_type": adata.obs["role"],
})
metadata_df.to_csv(PROJECT_ROOT / "data/metadata.txt", sep="\t", index=False)

# Map MGI -> HGNC gene symbols
symbol_map = (
    interaction_df[["mgi_symbol", "hgnc_symbol"]]
    .drop_duplicates()
    .set_index("mgi_symbol")["hgnc_symbol"]
    .to_dict()
)
adata.var["hgnc_symbol"] = adata.var_names.map(symbol_map)
adata = adata[:, adata.var["hgnc_symbol"].notna()].copy()
adata.var_names = adata.var["hgnc_symbol"]

# Export expression matrix
adata.write(PROJECT_ROOT / "data/mouse_e18.5d_counts.h5ad")

# Run CellPhoneDB
cpdb_statistical_analysis_method.call(
    cpdb_file_path=str(NGS_ROOT / "resources/cellphonedb/v5.0.0/cellphonedb.zip"),
    meta_file_path=str(PROJECT_ROOT / "data/metadata.txt"),
    counts_file_path=str(PROJECT_ROOT / "data/mouse_e18.5d_counts.h5ad"),
    counts_data="hgnc_symbol",
    output_path=str(PROJECT_ROOT / "output/cellphonedb"),
    score_interactions=True,
    threads=16,
)
