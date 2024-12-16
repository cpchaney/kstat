library(tidyverse)
library(Matrix)
library(ggpubr)
library(rjson)
library(here)
library(viridis)
library(pheatmap)
library(plot3D)

# Utility to load JSON configuration
load_config <- function(age) {
  here("config", str_c("stat_", age, "_16x16.json")) %>%
    fromJSON(file = .)
}

# Read and preprocess single-cell experiment data
load_sce <- function(age) {
  here("../foxd1-cre_myc_mycn/data", str_c("mouse_", age, "_foxd1-cre_myc_mycn_sce.rds")) %>%
    readRDS() %>%
    subset(, condition %in% c("wild type"))
}

# Generate color bar for visualization
generate_color_bar <- function(palette, file_path) {
  tibble(x = c(0, 1), y = c(0, 1)) %>%
    ggplot(aes(x = x, y = y, colour = y)) +
    geom_point() +
    scale_color_gradientn(name = "", colours = palette, limits = c(0, 1)) +
    theme_minimal() %>%
    ggpubr::as_ggplot(ggpubr::get_legend(.)) %>%
    ggsave(filename = file_path, device = "tiff", width = 2.125, height = 2.125, units = "in")
}

# Load spatial signal data
load_spatial_data <- function(age) {
  list(
    probabilities = here("data", str_c("mouse_", age, "_stat_cell_probabilities.mtx")) %>% readMM(),
    unique_bins = here("data", str_c("mouse_", age, "_stat_unique_bins.mtx")) %>% readMM(),
    unique_bins_key = here("data", str_c("mouse_", age, "_stat_unique_bins_key.csv")) %>% read_csv() %>%
      mutate(across(everything(), ~ . + 1))
  )
}

# Read and filter landmark data
load_landmarks <- function(age, expression_matrix) {
  here("data", str_c("cartana_", age, "_landmarks_meta.csv")) %>%
    read_csv() %>%
    pull(gene) %>%
    intersect(row.names(expression_matrix))
}

# Compute coarse probabilities for tiles
compute_coarse_probabilities <- function(tiles, probabilities, unique_bins_key) {
  do.call(cbind, lapply(tiles, function(tile) {
    rowMeans(as.matrix(probabilities[, unique_bins_key$index][, as.numeric(t(tile))]))
  }))
}

# Normalize cluster expression matrix
normalize_cluster_expression <- function(expression_matrix) {
  expression_matrix %>%
    as.data.frame() %>%
    rowwise() %>%
    mutate(max_value = max(c_across())) %>%
    ungroup() %>%
    mutate(across(everything(), ~ . / max_value)) %>%
    select(-max_value) %>%
    as_tibble()
}

# Generate heatmap
generate_heatmap <- function(matrix, color_palette, file_path, height = 8, width = 8) {
  pheatmap(
    matrix,
    color = color_palette,
    breaks = seq(0, 1, by = 0.01),
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    cellheight = height,
    cellwidth = width,
    fontsize_row = 8,
    legend = FALSE,
    annotation_legend = TRUE,
    annotation_names_col = FALSE
  ) %>%
    ggsave(filename = file_path, device = "tiff", width = width, height = height, units = "in")
}

# Save spatial signals as heatmap
save_spatial_signal <- function(spatial_signal, file_path, height, width, palette, config) {
  blurred_signal <- gaussianBlur(spatial_signal[unique_bins_key$index],
                                 height = config$targetHeight,
                                 width = config$targetWidth)

  spatial_matrix <- matrix(blurred_signal, nrow = config$targetHeight, byrow = TRUE)
  coarse_matrix <- matrix(sapply(split_matrix_into_tiles(spatial_matrix, 16), mean),
                          byrow = TRUE,
                          ncol = ceiling(config$targetWidth / 16)) %>%
    normalize_matrix()

  spatial_matrix %>%
    getSignalCoordinates() %>%
    ggplot(aes(x = X, y = Y)) +
    geom_raster(aes(fill = magnitude), interpolate = TRUE) +
    scale_fill_viridis(name = "Spatial Signal", option = "B") +
    coord_fixed(xlim = c(1, width), ylim = c(1, height), expand = FALSE) +
    theme_dark() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    ) +
    ggsave(filename = file_path, device = "tiff", width = 8.5, height = 11, units = "in")
}

# Main script workflow
main <- function(age = "e18.5d") {
  config <- load_config(age)
  sce <- load_sce(age)

  expression_matrix <- sce %>%
    assay("imputed") %>%
    .[rowSums(.) > 0, ]

  spatial_data <- load_spatial_data(age)
  landmarks <- load_landmarks(age, expression_matrix)

  cluster_probabilities <- compute_cluster_probabilities(sce, spatial_data$probabilities)
  coarse_cluster_probabilities <- compute_coarse_probabilities(
    tiles = split_matrix_into_tiles(config$template, 16),
    probabilities = cluster_probabilities,
    unique_bins_key = spatial_data$unique_bins_key
  )

  cluster_expression <- aggregate_cluster_expression(sce, expression_matrix, landmarks)
  normalized_expression <- normalize_cluster_expression(cluster_expression)

  generate_heatmap(
    normalized_expression,
    color_palette = magma(100),
    file_path = here("image/stat_cluster_expression.tiff")
  )

  save_spatial_signal(
    spatial_signal = getSpatialSignal(
      as.numeric(normalized_expression[landmarks == "Smoc2", ]),
      coarse_cluster_probabilities,
      NULL
    ),
    file_path = here("image/stat_expected_expression_smoc2.tiff"),
    height = config$targetHeight %/% 16 + 1,
    width = config$targetWidth %/% 16 + 1,
    palette = magma(100),
    config = config
  )
}

# Run the main function
main()
