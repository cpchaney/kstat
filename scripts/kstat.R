# Load required libraries
library(tidyverse)
library(rjson)
library(here)
library(patchwork)
library(purrr)

# Load configuration
load_config <- function(stage) {
  here("config", str_c("mouse_", stage, ".json")) %>% fromJSON(file = .)
}

# Load expression data
load_expression_data <- function(file_path) {
  readRDS(file = file_path) %>% assay("imputed")
}

# Load matrices and data
load_landmarks_and_bins <- function(stage, config) {
  list(
    landmarks_count_matrix = Matrix::readMM(here("data", str_c("mouse_", stage, "_landmarks_matrix.mtx"))),
    cells_bins_probabilities = Matrix::readMM(here("data", str_c(config$experiment_name, "_cells_bins_probabilities.mtx"))),
    unique_bins_key = read.csv(here("data", str_c(config$experiment_name, "_unique_bins_key.csv"))) %>%
      mutate(index = index + 1) %>% pull(index)
  )
}

# Generate a spatial plot for a single landmark
generate_landmark_plot <- function(gene_of_interest, landmark_signal, width, height, config) {
  signal_table <- getSignalCoordinates(
    gene_of_interest, landmark_signal, height = height, width = width
  )

  getSpatialPlot(
    signal_table,
    width = width,
    height = height,
    point_size = 0.2,
    channel_color_scale = setNames(c("red"), c(gene_of_interest)),
    background_file = here("image", str_c(config$experiment_name, "_background_scaled.jpg")),
    plot_legend = FALSE
  ) +
    annotate(
      "text",
      x = width - 24, y = 24,
      label = gene_of_interest,
      hjust = 1, vjust = 0,
      color = "white", size = 2, fontface = "bold"
    )
}

# Generate combined spatial plot for multiple landmarks
generate_combined_landmark_plot <- function(landmarks, landmarks_matrix, width, height, config) {
  plots <- sapply(landmarks, function(gene) {
    signal <- landmarks_matrix[gene, ]
    generate_landmark_plot(gene, signal, width, height, config)
  }, simplify = FALSE, USE.NAMES = TRUE)

  Reduce(`+`, plots) +
    plot_layout(ncol = 15) &
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
}

# Save combined landmark plot
save_combined_plot <- function(plot, filename, composite_width, composite_height, scale_factor = 8.5) {
  ggsave(
    filename = filename,
    plot = plot,
    device = "tiff",
    width = scale_factor * composite_width,
    height = scale_factor * composite_height,
    dpi = 300,
    bg = "transparent"
  )
}

# Generate thresholded signal plots for multiple genes
generate_thresholded_signal_plots <- function(genes, expression_matrix, probabilities, unique_bins_key, width, height, config) {
  signals <- do.call(rbind, lapply(genes, function(gene) {
    getSpatialSignal(expression_matrix[gene, ], probabilities, mapped.bins = NULL, log = FALSE)
  }))
  row.names(signals) <- genes
  signals <- signals[, unique_bins_key]

  plotThresholdedSignals(
    signals,
    height = height,
    width = width,
    percentiles = rep(95, length(genes)),
    denoise = TRUE,
    point_size = 0.2,
    col = c("red", "yellow", "white")[1:length(genes)],
    plot_legend = FALSE
  )
}

# Main function to generate visualizations and save results
main <- function(stage = "e18.5d") {
  # Load configuration
  config <- load_config(stage)

  # Load data
  expression_matrix <- load_expression_data("../kidney_development/data/mouse_e18.5d_snuc_sce.rds")
  data <- load_landmarks_and_bins(stage, config)
  width <- config$scaled_width
  height <- config$scaled_height

  # Plot measured expression for landmarks
  ordered_landmarks <- read.csv("data/mouse_e18.5d_cartana_genes_ordered.csv")$gene
  ordered_landmarks <- c(ordered_landmarks, setdiff(config$landmarks, ordered_landmarks))

  combined_landmark_plot <- generate_combined_landmark_plot(
    landmarks = ordered_landmarks,
    landmarks_matrix = data$landmarks_count_matrix,
    width = width,
    height = height,
    config = config
  )

  # Save combined plot
  composite_width <- 15 * (width / 300)
  composite_height <- 6 * (height / 300)
  save_combined_plot(
    combined_landmark_plot,
    here("image", str_c(config$experiment_name, "_combined_measured_landmarks.tiff")),
    composite_width,
    composite_height
  )

  # Plot imputed expression for a gene
  gene_of_interest <- "Gucy1a1"
  spatial_signal <- getSpatialSignal(
    as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
    data$cells_bins_probabilities,
    mapped.bins = NULL,
    log = FALSE
  )

  p <- plotExpectedSignal(
    spatial_signal[data$unique_bins_key],
    gene_of_interest,
    height,
    width,
    denoise = TRUE,
    plot_legend = FALSE
  )

  ggsave(
    filename = here("image", str_c(config$experiment_name, "_", tolower(gene_of_interest), "_expected.tiff")),
    plot = p,
    width = width,
    height = height,
    units = "px",
    dpi = 300
  )

  # Thresholded signal plot for multiple genes
  genes_of_interest <- c("Gucy1a2", "Clca3a1", "Cldn11")
  thresholded_signals_plot <- generate_thresholded_signal_plots(
    genes_of_interest,
    expression_matrix,
    data$cells_bins_probabilities,
    data$unique_bins_key,
    width,
    height,
    config
  )

  ggsave(
    filename = here("image", str_c("mouse_", stage, "_kidney_combined_inferred_", paste(tolower(genes_of_interest), collapse = "_"), ".tiff")),
    plot = thresholded_signals_plot,
    device = "tiff",
    height = height,
    width = width,
    units = "px",
    dpi = 300
  )
}
