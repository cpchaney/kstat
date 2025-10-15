load_interaction_data <- function(interaction_id, config, cell_types) {
  hdf5_path <- file.path("data", sprintf(
          "mouse_%s_kidney_%s_transport_distributions.h5",
          config$stage, paste0(cell_types, collapse = "_")
          ))
  print(hdf5_path)

  read_tensor <- function(path) as.numeric(h5read(hdf5_path, path))
  read_indices <- function(path) h5read(hdf5_path, path)

  list(
    ligand   = read_tensor(paste0("/", interaction_id, "/ligand")),
    receptor = read_tensor(paste0("/", interaction_id, "/receptor")),
    delivered = read_tensor(paste0("/", interaction_id, "/delivered")),
    received  = read_tensor(paste0("/", interaction_id, "/received")),
    delivered_idx = read_indices(paste0("/", interaction_id, "/delivered_indices")),
    received_idx  = read_indices(paste0("/", interaction_id, "/received_indices"))
  )
}
