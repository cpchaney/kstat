#' Initialize the R analysis environment.
#'
#' Loads required R packages (with suppressed startup messages) and sets
#' the active Python conda environment via `reticulate`.
#'
#' @param env_name Character string naming the conda environment to activate (default = "statlas").
#'
#' @return None. Side effects: loads libraries and sets conda environment.
#'
#' @examples
#' initialize_environment("my_env")
initialize_environment <- function(env_name = "statlas") {
  suppressPackageStartupMessages({
    # Load core packages for analysis, visualization, and data handling
    library(tidyverse)  # Includes dplyr, ggplot2, etc.
    library(rjson)      # For JSON config files
    library(here)       # Simplifies file paths
    library(rhdf5)      # HDF5 file access
    library(cowplot)    # ggplot2 extensions for layout
    library(patchwork)  # ggplot2 composition
    library(ggpubr)     # Publication-ready plots
    library(Matrix)     # Sparse matrix support
    library(ragg)       # High-quality graphics device
  })

  # Activate specified Python conda environment
  reticulate::use_condaenv(env_name, required = TRUE)

  # Confirmation message
  message("Environment initialized: ", env_name)
}
