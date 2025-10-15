#!/usr/bin/env Rscript

# ---------------------- #
# Setup & Configuration #
# ---------------------- #

library(tidyverse)
library(rjson)
library(here)
library(Matrix)
library(patchwork)

# Source local utilities
source(here("R/utils.R"))
source(here("R/visualization.R"))

# Load config
stage <- "e18.5d"
config_path <- here("config", paste0("mouse_", stage, ".json"))
config <- fromJSON(file = config_path)

# Define shared variables
width <- config$scaled_width
height <- config$scaled_height
landmarks <- config$landmarks

background_file <- here("image", paste0(config$experiment_name, "_background_scaled.jpg"))
cell_bin_file <- here("data", paste0(config$experiment_name, "_cells_bins_probabilities.mtx"))
bin_key_file <- here("data", paste0(config$experiment_name, "_unique_bins_key.csv"))

# ------------------------- #
# Load Data                #
# ------------------------- #

# Landmark matrix
landmark_path <- here("data", paste0("mouse_", stage, "_landmarks_matrix.mtx"))
landmark_counts <- Matrix::readMM(landmark_path)
rownames(landmark_counts) <- landmarks

# Cell-bin probabilities
cells_bins_probabilities <- Matrix::readMM(cell_bin_file)

# Bin mapping (convert 0-based to 1-based index)
unique_bins_key <- read.csv(bin_key_file)$index + 1

# ------------------------- #
# Function: Save Plot      #
# ------------------------- #

save_expression_plot <- function(plot, identifier, suffix = "expression") {
  out_file <- here("image", sprintf(
    "%s_%s_%s.tiff",
    config$experiment_name,
    stringr::str_to_lower(identifier),
    suffix
  ))

  ggplot2::ggsave(
    filename = out_file,
    plot = plot + theme(text = element_text(family = "Arial")),
    device = "tiff",
    height = height,
    width = width,
    units = "px",
    dpi = 300,
    compression = "lzw"
  )
  message("? Saved plot:", out_file)
}

# ------------------------- #
# 1. Plot Measured Landmark #
# ------------------------- #

gene_of_interest <- "Mgp"
signal_vector <- landmark_counts[gene_of_interest, ]

signal_table <- getSignalCoordinates(
  channel_name = gene_of_interest,
  v = signal_vector,
  height = height,
  width = width
)

measured_plot <- getSpatialPlot(
  signal_table,
  width = width,
  height = height,
  point_size = 0.5,
  channel_color_scale = setNames("red", gene_of_interest),
  background_file = background_file,
  plot_legend = FALSE
)

print(measured_plot)
save_expression_plot(measured_plot, gene_of_interest, "measured_expression")

# ------------------------- #
# 2. Plot Inferred Gene     #
# ------------------------- #

gene_of_interest <- "Cspg4"

spatial_signal <- getSpatialSignal(
  signal = as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
  cell_probabilities = cells_bins_probabilities,
  mapped.bins = NULL,
  log = FALSE
)

inferred_plot <- plotThresholdedSignal(
  spatial_signal = spatial_signal[unique_bins_key],
  signal_identifier = gene_of_interest,
  height = height,
  width = width,
  point_size = 0.5,
  denoise = TRUE,
  percentile = 96,
  background_file = background_file,
  plot_legend = FALSE
)

print(inferred_plot)
save_expression_plot(inferred_plot, gene_of_interest, "inferred_expression")

# ------------------------- #
# 3. Plot Regulon Activity  #
# ------------------------- #

regulon_of_interest <- "Lef1"

regulons_auc_matrix <- t(metadata(sce)[["regulons_auc_matrix"]])[, colnames(sce)]

spatial_signal <- getSpatialSignal(
  signal = as.numeric(regulons_auc_matrix[regulon_of_interest, , drop = FALSE]),
  cell_probabilities = cells_bins_probabilities,
  mapped.bins = NULL,
  log = FALSE
)

regulon_plot <- plotThresholdedSignal(
  spatial_signal = spatial_signal[unique_bins_key],
  signal_identifier = regulon_of_interest,
  height = height,
  width = width,
  point_size = 0.5,
  denoise = TRUE,
  percentile = 96,
  background_file = background_file,
  plot_legend = FALSE
)

print(regulon_plot)
save_expression_plot(regulon_plot, regulon_of_interest, "inferred_activity")

# ----------------------------- #
# 4. Plot Gene Set Enrichment  #
# ----------------------------- #

gene_sets_auc_matrix <- metadata(sce)[["gene_sets_auc_matrix"]][, colnames(sce)]
gene_set_of_interest <- "HALLMARK_G2M_CHECKPOINT"

spatial_signal <- getSpatialSignal(
  signal = as.numeric(gene_sets_auc_matrix[gene_set_of_interest, , drop = FALSE]),
  cell_probabilities = cells_bins_probabilities,
  mapped.bins = NULL,
  log = FALSE
)

gene_set_plot <- plotThresholdedSignal(
  spatial_signal = spatial_signal[unique_bins_key],
  signal_identifier = gene_set_of_interest,
  height = height,
  width = width,
  point_size = 0.5,
  denoise = TRUE,
  percentile = 98,
  background_file = background_file,
  plot_legend = FALSE
)

print(gene_set_plot)

# Safe filename (first 3 letters of each component)
gene_set_label <- gene_set_of_interest |>
  strsplit("_") |>
  unlist() |>
  substr(1, 3) |>
  paste0(collapse = "_") |>
  tolower()

save_expression_plot(gene_set_plot, gene_set_label, "inferred_activity")
