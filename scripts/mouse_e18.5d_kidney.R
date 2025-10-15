library(tidyverse)
library(rjson)
library(here)
library(patchwork)
library(purrr)
library(rhdf5)

# Source utility scripts
source(here("../tools/R/sc-rna-seq_util.R"))
source(here("R/utils.R"))
source(here("R/visualization.R"))

use_condaenv("statlas", required = TRUE)

# Define age and load configuration
stage <- "e18.5d"
config <- here("config", str_c("mouse_", stage, ".json")) %>%
  fromJSON(file = .)

sce <- readRDS(file = "../kidney_development/mouse_e18.5d/data/mouse_e18.5d_sn_sce.rds")
expression_matrix <- assay(sce, "imputed")

width <- config$scaled_width
height <- config$scaled_height

# Load the data
landmarks_count_matrix <- Matrix::readMM(paste0("data/mouse_", stage, "_landmarks_matrix.mtx"))
cells_bins_probabilities <- Matrix::readMM(paste0(
  "data/",
  config$experiment_name,
  "_cells_bins_probabilities.mtx"
))
unique_bins_key <- read.csv(file = paste0(
  "data/",
  config$experiment_name,
  "_unique_bins_key.csv"
))$index + 1
common_landmarks <- config$common_landmarks

# Prepare data
landmarks <- config$landmarks
row.names(landmarks_count_matrix) <- landmarks

ordered_landmarks <- landmarks[order(landmarks)]

# Plot measured expression
landmark_plots <- sapply(ordered_landmarks, function(gene_of_interest) {
        landmark_signal <- landmarks_count_matrix[gene_of_interest, ]

        landmark_signal_table <- getSignalCoordinates(gene_of_interest,
                landmark_signal,
                height = height,
                width = width
        )

        return(
                getSpatialPlot(
                        landmark_signal_table,
                        width = width,
                        height = height,
                        point_size = 0.5,
                        channel_color_scale = setNames(c("red"), c(gene_of_interest)),
                        background_file = paste0(
                                "image/",
                                config$experiment_name,
                                "_background_scaled.jpg"
                        ),
                        plot_legend = FALSE
                ) +
                        annotate("text",
                                x = width - 24,
                                y = 0 + 24,
                                label = gene_of_interest,
                                hjust = 1,
                                vjust = 0,
                                color = "white",
                                size = 2,
                                fontface = "bold"
                        )
        )
}, simplify = FALSE, USE.NAMES = TRUE)

combined_landmark_plot <- Reduce(`+`, landmark_plots) +
  plot_layout(ncol = 8) &
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    element_text(family = "Arial"),
    panel.background = element_rect(fill = "black", color = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 0,
      # Top margin
      r = 0,
      # Right margin
      b = 0,
      # Bottom margin
      l = 0
    )
  )

composite_width <- 8 * (width / 300)
composite_height <- 12 * (height / 300)
scale_factor <- 8.5 / composite_width

ggsave(
  filename = paste0(
    "image/",
    config$experiment_name,
    "_combined_measured_landmarks.tiff"
  ),
  # Output file name
  plot = combined_landmark_plot,
  # Plot object
  device = "tiff",
  # Save as TIFF
  width = scale_factor * composite_width,
  # Total width in inches
  height = scale_factor * composite_height,
  # Total height in inches
  dpi = 300, # Resolution in DPI
  bg = "transparent"
)

# Plot measuured expression of individual landmarks
gene_of_interest <- "Mgp"

landmark_signal <- landmarks_count_matrix[gene_of_interest, ]

landmark_signal_table <- getSignalCoordinates(gene_of_interest,
  landmark_signal,
  height = height,
  width = width
)

p <- getSpatialPlot(
        landmark_signal_table,
        width = width,
        height = height,
        point_size = 0.5,
        channel_color_scale = setNames(c("red"), c(gene_of_interest)),
        background_file = paste0(
                "image/",
                config$experiment_name,
                "_background_scaled.jpg"
        ),
        plot_legend = F
)

ggsave(
  filename = paste0(
    "image/",
    config$experiment_name,
    "_",
    stringr::str_to_lower(gene_of_interest),
    "_measured.tiff"
  ),
  plot = p,
  width = width,
  height = height,
  units = "px",
  bg = "transparent",
  dpi = 300
)

# Plot imputed expression
gene_of_interest <- "Cspg4"
spatial_signal <- getSpatialSignal(as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
  cells_bins_probabilities,
  mapped.bins = NULL,
  log = F
)

p <- plotExpectedSignal(spatial_signal[unique_bins_key],
  gene_of_interest,
  height,
  width,
  denoise = F,
  plot_legend = FALSE
)
p
ggsave(
  filename = paste0(
    "image/",
    config$experiment_name,
    "_",
    stringr::str_to_lower(gene_of_interest),
    "_expected.tiff"
  ),
  plot = p,
  width = width,
  height = height,
  units = "px",
  bg = "transparent",
  dpi = 300
)

threshold <- 95
p <- plotThresholdedSignal(
  spatial_signal[unique_bins_key],
  point_size = 0.5,
  gene_of_interest,
  height = height,
  width = width,
  background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
  percentile = threshold,
  plot_legend = FALSE
)
p
ggsave(
  filename = paste0(
    "image/",
    config$experiment_name,
    "_",
    stringr::str_to_lower(gene_of_interest),
    "_inferred_",
    threshold,
    ".tif"
  ),
  plot = p,
  width = width,
  height = height,
  units = "px",
  bg = "transparent",
  dpi = 300
)


common_nonzero_landmarks <- common_landmarks[rowSums(expression_matrix[common_landmarks, ]) != 0]

inferred_landmark_plots <- sapply(common_nonzero_landmarks, function(gene_of_interest) {
  spatial_signal <- getSpatialSignal(as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
    cells_bins_probabilities,
    mapped.bins = NULL,
    log = F
  )
  p <- plotThresholdedSignal(
    spatial_signal[unique_bins_key],
    gene_of_interest,
    height = height,
    width = width,
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    percentile = 99,
    plot_legend = FALSE,
    point_size = 0.2
  )

  return(
    p +
      annotate("text",
        x = -Inf,
        y = Inf,
        label = gene_of_interest,
        hjust = -0.1,
        vjust = 1.5,
        size = 2,
        fontface = "bold",
        color = "white"
      ) +
      theme(element_text(family = "Arial"))
  )
}, simplify = FALSE, USE.NAMES = TRUE)

combined_inferred_landmarks_plot <- Reduce(`+`, inferred_landmark_plots) +
  plot_layout(ncol = 8) &
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    plot.background = element_rect(fill = "black", color = NA), # Entire plot background
    panel.background = element_rect(fill = "black", color = NA)
  )

ggsave(
        filename = paste0(
                "image/",
                config$experiment_name,
                "_combined_inferred_landmarks.tiff"
        ),
        combined_inferred_landmarks_plot,
        device = "tiff",
        height = config$scaled_height,
        width = config$scaled_width,
        units = "px"
)

# Figure 2 generation
genes_of_interest <- c("Podxl", "Smoc2", "Slc27a2", "Egfl7")
gene_panels <- sapply(genes_of_interest, function(gene_of_interest) {
  panel_plots <- list()
  landmark_signal <- landmarks_count_matrix[gene_of_interest, ]

  landmark_signal_table <- getSignalCoordinates(gene_of_interest,
    landmark_signal,
    height = height,
    width = width
  )

  p <- getSpatialPlot(
    landmark_signal_table,
    width = width,
    height = height,
    point_size = 0.1,
    channel_color_scale = setNames(c("red"), c(gene_of_interest)),
    background_file = paste0(
      "image/",
      config$experiment_name,
      "_background_scaled.jpg"
    ),
    plot_legend = F
  )
  panel_plots <- c(panel_plots, list(p))

  spatial_signal <- getSpatialSignal(as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
    cells_bins_probabilities,
    mapped.bins = NULL,
    log = F
  )

  p <- plotExpectedSignal(spatial_signal[unique_bins_key],
    gene_of_interest,
    height,
    width,
    denoise = T,
    plot_legend = FALSE
  )
  panel_plots <- c(panel_plots, list(p))

  threshold_panels <- lapply(c(90, 95, 99), function(threshold) {
    p <- plotThresholdedSignal(
      spatial_signal[unique_bins_key],
      point_size = 0.,
      gene_of_interest,
      height = height,
      width = width,
      background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
      percentile = threshold,
      plot_legend = FALSE
    )
    return(p)
  })
  panel_plots <- c(panel_plots, threshold_panels)

  combined_gene_panel_plot <- Reduce(`+`, panel_plots) +
    plot_layout(ncol = 1, widths = 1) &
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      plot.background = element_rect(fill = "black", color = NA), # Entire plot background
      panel.background = element_rect(fill = "black", color = NA)
    )

  ggsave(
    filename = paste0(
      "image/",
      config$experiment_name,
      "_",
      stringr::str_to_lower(gene_of_interest),
      "_paneled.tiff"
    ),
    plot = combined_gene_panel_plot,
    height = 3.25,
    width = 0.70,
    units = "in",
    bg = "transparent",
    dpi = 300
  )

  return(combined_gene_panel_plot)
}, simplify = FALSE, USE.NAMES = TRUE)

# Figure 3 generation
genes_of_interest <- c("Pou3f3", "Ass1", "Cldn11", "Apcdd1")

inferred_gene_plots <- sapply(genes_of_interest, function(gene_of_interest) {
  spatial_signal <- getSpatialSignal(as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
    cells_bins_probabilities,
    mapped.bins = NULL,
    log = F
  )
  p <- plotThresholdedSignal(
    spatial_signal[unique_bins_key],
    gene_of_interest,
    height = height,
    width = width,
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    percentile = 99,
    plot_legend = FALSE,
    point_size = 0.2
  ) +
    annotate("text",
      x = height - 24,
      y = 0 + 24,
      label = gene_of_interest,
      hjust = 1,
      vjust = 0,
      color = "white",
      size = 4,
      fontface = "bold"
    )

  return(p)
}, simplify = FALSE, USE.NAMES = TRUE)

combined_landmark_plot <- Reduce(`+`, inferred_gene_plots) +
  plot_layout(ncol = 4) &
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    element_text(family = "Arial"),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "black", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "black", color = "black"),
    plot.margin = margin(
      t = 0,
      # Top margin
      r = 0,
      # Right margin
      b = 0,
      # Bottom margin
      l = 0
    )
  )

composite_width <- 4 * (height / 300)
composite_height <- 1 * (width / 300)
scale_factor <- 8.5 / composite_width

ggsave(
  filename = paste0(
    "image/",
    config$experiment_name,
    "_combined_inferred_",
    paste0(stringr::str_to_lower(genes_of_interest), collapse = "_"),
    ".tiff"
  ),
  # Output file name
  plot = combined_landmark_plot,
  # Plot object
  device = "tiff",
  # Save as TIFF
  width = scale_factor * composite_width,
  # Total width in inches
  height = scale_factor * composite_height,
  units = "in",
  # Total height in inches
  dpi = 300, # Resolution in DPI
  bg = "transparent"
)

cluster_means <- cbind(
  as.data.frame(as.matrix(t(expression_matrix[common_landmarks, ]))),
  as.data.frame(colData(sce)[, "cluster", drop = FALSE])
) %>%
  group_by(cluster) %>%
  summarize_all(mean) %>%
  column_to_rownames(var = "cluster")

landmark_variances <- setNames(colVars(as.matrix(cluster_means)), common_landmarks)
landmark_variances[order(landmark_variances)]

landmark_means <- setNames(colMeans(as.matrix(cluster_means)), common_landmarks)
landmark_means[order(landmark_means, decreasing = TRUE)]

landmark_lower_quartiles <- setNames(rowQuantiles(expression_matrix[common_landmarks, ], probs = 0.25), common_landmarks)
landmark_lower_quartiles[order(landmark_lower_quartiles, decreasing = TRUE)]

landmark_table <- data.frame(
  landmark = common_landmarks,
  mean = landmark_means,
  cluster_variance = landmark_variances
)

ggplot(landmark_table, aes(x = mean, y = cluster_variance)) +
  geom_point(aes(color = ifelse(landmark %in% c("H19", "Fxyd2", "Acta2", "Ldhb", "Cdkn1c"), "red", "black")), size = 2) +
  geom_text(
    aes(label = ifelse(landmark %in% c("H19", "Fxyd2", "Acta2", "Ldhb", "Cdkn1c"), landmark, "")), # Add labels for points meeting the condition
    vjust = -1 # Adjust vertical position of labels
  ) +
  scale_color_identity() +
  theme_minimal()


landmark_quantiles <- data.frame(rowQuantiles(expression_matrix[common_landmarks, ], probs = c(0.25, 0.5, 0.75)))
colnames(landmark_quantiles) <- c("lower", "median", "upper")
landmark_quantiles <- rownames_to_column(landmark_quantiles, var = "landmark")

ggplot(landmark_quantiles, aes(x = median, y = landmark)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
  labs(
    title = "Interquartile range for landmark expression",
    x = "Expression",
    y = "Landmark"
  ) +
  theme_bw() #+
theme(
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  axis.line = element_line(colour = "black"),
  legend.position = "none",
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 20, face = "bold")
)


################################
## Visualize regulon activity ##
################################
regulons_auc_matrix <- t(metadata(sce)[["regulons_auc_matrix"]])[, colnames(sce)]
regulon_of_interest <- "Lef1"

spatial_signal <- getSpatialSignal(as.numeric(regulons_auc_matrix[regulon_of_interest, , drop = FALSE]),
  cells_bins_probabilities,
  mapped.bins = NULL,
  log = F
)

plotExpectedSignal(spatial_signal[unique_bins_key],
  regulon_of_interest,
  height,
  width,
  denoise = T
)


p <- plotThresholdedSignal(
  spatial_signal[unique_bins_key],
  point_size = 0.5,
  regulon_of_interest,
  height = height,
  width = width,
  background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
  percentile = 96,
  plot_legend = FALSE
)
p

ggsave(
  filename = paste0(
    "image/",
    config$experiment_name,
    "_",
    stringr::str_to_lower(regulon_of_interest),
    "_inferred_acitivity.tiff"
  ),
  plot = p + theme(element_text(family = "Arial")),
  device = "tiff",
  height = 11,
  width = 8.5,
  units = "in"
)

################################
## Visualize gene set enrichment ##
################################
gene_sets_auc_matrix <- metadata(sce)[["gene_sets_auc_matrix"]][, colnames(sce)]
gene_set_of_interest <- "HALLMARK_G2M_CHECKPOINT"

spatial_signal <- getSpatialSignal(as.numeric(gene_sets_auc_matrix[gene_set_of_interest, , drop = FALSE]),
  cells_bins_probabilities,
  mapped.bins = NULL,
  log = F
)

plotExpectedSignal(spatial_signal[unique_bins_key],
  gene_set_of_interest,
  height,
  width,
  denoise = T
)


p <- plotThresholdedSignal(
  spatial_signal[unique_bins_key],
  point_size = 0.5,
  gene_set_of_interest,
  height = height,
  width = width,
  background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
  percentile = 98,
  plot_legend = FALSE
)
p

ggsave(
  filename = paste0(
    "image/",
    config$experiment_name,
    "_",
    stringr::str_to_lower(paste(substr(strsplit(gene_set_of_interest, "_")[[1]], 1, 3), collapse = "_")),
    "_inferred_acitivity.tiff"
  ),
  plot = p + theme(element_text(family = "Arial")),
  device = "tiff",
  height = config$scaled_height,
  width = config$scaled_width,
  units = "px"
)


