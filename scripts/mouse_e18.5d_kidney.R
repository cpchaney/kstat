# ---- Setup & Configuration ----
library(tidyverse)
library(rjson)
library(here)
library(patchwork)
library(purrr)
library(rhdf5)

# Load utility scripts
source(here("../tools/R/sc-rna-seq_util.R"))
source(here("R/utils.R"))
source(here("R/visualization.R"))

use_condaenv("kstat", required = TRUE)

stage <- "e18.5d"
config <- fromJSON(file = here("config", str_c("mouse_", stage, ".json")))

# ---- Data Loading ----
sce <- readRDS(file = "../kidney_development/mouse_e18.5d/data/mouse_e18.5d_sn_sce.rds")
expression_matrix <- assay(sce, "imputed")
width <- config$scaled_width
height <- config$scaled_height

landmarks_count_matrix <- Matrix::readMM(paste0("data/mouse_", stage, "_landmarks_matrix.mtx"))
cells_bins_probabilities <- Matrix::readMM(paste0("data/", config$experiment_name, "_cells_bins_probabilities.mtx"))
unique_bins_key <- read.csv(paste0("data/", config$experiment_name, "_unique_bins_key.csv"))$index + 1

landmarks <- config$landmarks
row.names(landmarks_count_matrix) <- landmarks
ordered_landmarks <- sort(landmarks)
common_landmarks <- config$common_landmarks

# ---- Utility Functions ----
plot_measured_landmark <- function(gene, width, height) {
  signal <- landmarks_count_matrix[gene, ]
  signal_table <- getSignalCoordinates(gene, signal, height, width)

  getSpatialPlot(
    signal_table,
    width = width,
    height = height,
    point_size = 0.5,
    channel_color_scale = setNames("red", gene),
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    plot_legend = FALSE
  ) +
    annotate("text", x = width - 24, y = 24, label = gene, hjust = 1, vjust = 0, color = "white", size = 2, fontface = "bold")
}

plot_inferred_expression <- function(gene, threshold = 99, point_size = 0.2) {
  signal <- getSpatialSignal(as.numeric(expression_matrix[gene, , drop = FALSE]),
                              cells_bins_probabilities, mapped.bins = NULL, log = FALSE)
  plotThresholdedSignal(
    signal[unique_bins_key], gene, height, width,
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    percentile = threshold, plot_legend = FALSE, point_size = point_size
  ) +
    annotate("text", x = -Inf, y = Inf, label = gene, hjust = -0.1, vjust = 1.5, size = 2, fontface = "bold", color = "white") +
    theme(element_text(family = "Arial"))
}

save_plot <- function(plot, filename, width, height, scale_factor = 1, dpi = 300, units = "px") {
  ggsave(
    filename = filename,
    plot = plot,
    width = width * scale_factor,
    height = height * scale_factor,
    dpi = dpi,
    bg = "transparent",
    device = "tiff",
    units = units
  )
}

# ---- Composite Measured Landmark Plot ----
landmark_plots <- map(ordered_landmarks, ~plot_measured_landmark(.x, width, height))
combined_landmark_plot <- wrap_plots(landmark_plots, ncol = 8) &
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "black"),
    panel.border = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    element_text(family = "Arial")
  )

save_plot(
  combined_landmark_plot,
  paste0("image/", config$experiment_name, "_combined_measured_landmarks.tiff"),
  width = 8 * (width / 300),
  height = 12 * (height / 300),
  scale_factor = 8.5 / (8 * (width / 300)),
  units = "in"
)

# ---- Individual Measured Plot Example ----
plot_measured_landmark("Mgp", width, height) %>%
  save_plot(paste0("image/", config$experiment_name, "_mgp_measured.tiff"), width, height)

# ---- Individual Inferred Plot Example ----
plot_inferred_expression("Cspg4", threshold = 95) %>%
  save_plot(paste0("image/", config$experiment_name, "_cspg4_inferred_95.tif"), width, height)

# ---- Composite Inferred Landmark Plot ----
common_nonzero_landmarks <- common_landmarks[rowSums(expression_matrix[common_landmarks, ]) != 0]
inferred_landmark_plots <- map(common_nonzero_landmarks, ~plot_inferred_expression(.x))
combined_inferred_landmarks_plot <- wrap_plots(inferred_landmark_plots, ncol = 8) &
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA)
  )

save_plot(
  combined_inferred_landmarks_plot,
  paste0("image/", config$experiment_name, "_combined_inferred_landmarks.tiff"),
  config$scaled_width,
  config$scaled_height
)

# ---- Figure 2: Landmark Gene Panels ----
genes_of_interest <- c("Podxl", "Smoc2", "Slc27a2", "Egfl7")

map(genes_of_interest, function(gene) {
  panel_plots <- list(
    plot_measured_landmark(gene, width, height),
    plotExpectedSignal(
      getSpatialSignal(as.numeric(expression_matrix[gene, , drop = FALSE]), cells_bins_probabilities),
      gene, height, width, denoise = TRUE, plot_legend = FALSE
    )
  )

  threshold_plots <- map(c(90, 95, 99), ~plot_inferred_expression(gene, threshold = .x, point_size = 0))
  combined_panel <- wrap_plots(c(panel_plots, threshold_plots), ncol = 1) &
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      plot.background = element_rect(fill = "black"),
      panel.background = element_rect(fill = "black")
    )

  save_plot(
    combined_panel,
    paste0("image/", config$experiment_name, "_", str_to_lower(gene), "_paneled.tiff"),
    width = 0.7,
    height = 3.25,
    units = "in"
  )

  return(combined_panel)
})

# ---- Figure 3: Another Gene Panel ----
genes_of_interest <- c("Pou3f3", "Ass1", "Cldn11", "Apcdd1")
inferred_gene_plots <- map(genes_of_interest, plot_inferred_expression)

combined_landmark_plot <- wrap_plots(inferred_gene_plots, ncol = 4) &
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black"),
    plot.margin = margin(0, 0, 0, 0)
  )

save_plot(
  combined_landmark_plot,
  paste0("image/", config$experiment_name, "_combined_inferred_", paste0(str_to_lower(genes_of_interest), collapse = "_"), ".tiff"),
  width = 4 * (height / 300),
  height = 1 * (width / 300),
  scale_factor = 8.5 / (4 * (height / 300)),
  units = "in"
)

# ---- Landmark Expression Summary ----
cluster_means <- expression_matrix[common_landmarks, ] %>%
  t() %>%
  as.data.frame() %>%
  bind_cols(cluster = colData(sce)[["cluster"]]) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean)) %>%
  column_to_rownames("cluster")

landmark_means <- colMeans(cluster_means)
landmark_vars <- colVars(as.matrix(cluster_means))
landmark_table <- tibble(
  landmark = names(landmark_means),
  mean = landmark_means,
  cluster_variance = landmark_vars
)

ggplot(landmark_table, aes(x = mean, y = cluster_variance)) +
  geom_point(aes(color = ifelse(landmark %in% c("H19", "Fxyd2", "Acta2", "Ldhb", "Cdkn1c"), "red", "black")), size = 2) +
  geom_text(aes(label = ifelse(landmark %in% c("H19", "Fxyd2", "Acta2", "Ldhb", "Cdkn1c"), landmark, "")), vjust = -1) +
  scale_color_identity() +
  theme_minimal()

# ---- Interquartile Ranges ----
iqr_stats <- rowQuantiles(expression_matrix[common_landmarks, ], probs = c(0.25, 0.5, 0.75))
iqr_df <- as.data.frame(iqr_stats)
colnames(iqr_df) <- c("lower", "median", "upper")
iqr_df$landmark <- rownames(expression_matrix[common_landmarks, ])

ggplot(iqr_df, aes(x = median, y = landmark)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
  labs(title = "Interquartile range for landmark expression", x = "Expression", y = "Landmark") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# ---- Regulon Activity ----
regulons_auc_matrix <- t(metadata(sce)$regulons_auc_matrix)[, colnames(sce)]
regulon_of_interest <- "Lef1"

spatial_signal <- getSpatialSignal(as.numeric(regulons_auc_matrix[regulon_of_interest, , drop = FALSE]),
                                   cells_bins_probabilities)

plotExpectedSignal(spatial_signal[unique_bins_key], regulon_of_interest, height, width, denoise = TRUE)

plotThresholdedSignal(
  spatial_signal[unique_bins_key], point_size = 0.5, regulon_of_interest,
  height, width,
  background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
  percentile = 96, plot_legend = FALSE
) %>%
  save_plot(
    paste0("image/", config$experiment_name, "_", str_to_lower(regulon_of_interest), "_inferred_acitivity.tiff"),
    width = 8.5, height = 11, units = "in"
  )

# ---- Gene Set Enrichment ----
gene_sets_auc_matrix <- metadata(sce)$gene_sets_auc_matrix[, colnames(sce)]
gene_set <- "HALLMARK_G2M_CHECKPOINT"

spatial_signal <- getSpatialSignal(as.numeric(gene_sets_auc_matrix[gene_set, , drop = FALSE]),
                                   cells_bins_probabilities)

plotExpectedSignal(spatial_signal[unique_bins_key], gene_set, height, width, denoise = TRUE)

plotThresholdedSignal(
  spatial_signal[unique_bins_key], point_size = 0.5, gene_set,
  height, width,
  background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
  percentile = 98, plot_legend = FALSE
) %>%
  save_plot(
    paste0("image/", config$experiment_name, "_", str_to_lower(str_c(str_sub(str_split(gene_set, "_")[[1]], 1, 3), collapse = "_")), "_inferred_acitivity.tiff"),
    width = config$scaled_width, height = config$scaled_height
  )
