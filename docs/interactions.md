## 1. Configuration and Setup

The ligand-receptor inference pipeline is controlled via a YAML config file.

Example:

```yaml
project_root: /path/to/kstat/
ngs_root: /path/to/raw_ngs/
sql_query_path: queries/select_ligand_receptor.sql
database: cellchat.sqlite
```
We load the config and define project paths:
```python
from pathlib import Path
import yaml

with open("../../config/config.yaml") as f:
    config = yaml.safe_load(f)

PROJECT_ROOT = Path(config["project_root"])
NGS_ROOT = Path(config["ngs_root"])
```

---

## 2. Load Ligand-Receptor Table

We use a SQL query (stored in a file) to fetch curated ligand-receptor pairs from a SQLite database.

```python
with open(sql_path, "r") as f:
    sql_query = f.read()

with sqlite3.connect(db_path) as conn:
    interaction_df = pd.read_sql_query(sql_query, conn)
```
This generates a table like:

| id_cp_interaction | mgi_symbol | role     |
|-------------------|------------|----------|
| CPI-XYZ123        | Fgf9       | ligand   |
| CPI-XYZ123        | Fgfr1      | receptor |

Saved to:
`data/ligand_receptor_components.csv`

---

## 3. Load and Preprocess AnnData

We load a single-cell dataset and filter to genes involved in ligand-receptor interactions.

```python
adata = sc.read_h5ad("...")
# Filter for relevant genes
mgi_genes = pd.unique(interaction_df["mgi_symbol"].values)
adata = adata[:, adata.var_names.isin(mgi_genes)].copy()
```
We set `.X` to the imputed layer:

```python
adata.layers["counts"] = adata.X.copy()
adata.X = adata.layers["imputed"]
```
---

## 4. Load Cell-Bin Probabilities

We load the inferred probability matrix and bin mapping:

```python
cells_bins_probs = mmread("...cells_bins_probabilities.mtx").tocsr()
cells_bins_probs = torch.tensor(cells_bins_probs.toarray(), dtype=torch.float32)

unique_bins_key = torch.tensor(
    pd.read_csv("...unique_bins_key.csv")["index"]
)
```
- `cells_bins_probs`: tensor of shape (cells × bins)
- `unique_bins_key`: indexes valid spatial bins

---

## 5. Run the Interaction Inference Pipeline

We run `run_interaction_pipeline()` from the `kstat.interactions` module.

```python
results_df, distributions = run_interaction_pipeline(
    adata=adata,
    interaction_df=interaction_df,
    cell_type_list=["NPC1", "NPC2", "NPC3"],
    cell_autonomous=False,
    cells_bins_prob=cells_bins_probs,
    unique_bins_key=unique_bins_key,
    height=512,
    width=512,
    device="cuda",
    p=1,
    blur=3.0,
    reach=20.0,
    mode="unbalanced",
    return_distributions=True
)
```
- The output `results_df` contains interaction scores and metadata.
- The `distributions` object stores raw delivery/reception tensors per interaction.

---

## 6. Save Outputs

Results include:

- A summary CSV of interactions and scores
- Serialized receptor/ligand distributions (Pickle and HDF5)

```python
results_df.to_csv("output/..._interactions.csv", index=False)

with open("data/..._receptor_saturations.pkl", "wb") as f:
    pickle.dump(distributions, f)

save_distributions_to_hdf5(
    distributions,
    "data/..._transport_distributions.h5"
)
```
