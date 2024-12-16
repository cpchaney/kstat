library(tidyverse)
library(rjson)
library(here)
library(patchwork)
library(purrr)

# Source utility scripts
source(here("../tools/R/sc-rna-seq_util.R"))
source(here("R/utils.R"))
source(here("R/visualization.R"))

# Define age and load configuration
stage <- "e18.5d"
config <- here("config", str_c("mouse_", stage, ".json")) %>%
  fromJSON(file = .)

sce <- readRDS(file = "../kidney_development/data/mouse_e18.5d_snuc_sce.rds")
expression_matrix <- assay(sce, "imputed")

width <- config$scaled_width
height <- config$scaled_height

# Plot measured expression
landmarks_count_matrix <- Matrix::readMM(paste0("data/mouse_", stage, "_landmarks_matrix.mtx"))

landmarks <- config$landmarks
row.names(landmarks_count_matrix) <- landmarks

landmark_plots <- sapply(landmarks, function(gene_of_interest) {
  landmark_signal <- landmarks_count_matrix[gene_of_interest, ]
  
  landmark_signal_table <- getSignalCoordinates(gene_of_interest,
                                                landmark_signal,
                                                height = height,
                                                width = width)
  
  return(
    getSpatialPlot(
      landmark_signal_table,
      width = width,
      height = height,
      point_size = 0.2,
      channel_color_scale = setNames(c("red"), c(gene_of_interest)),
      background_file = paste0("image/",
                               config$experiment_name, 
                               "_background_scaled.jpg"),
      plot_legend = F
    ) + 
      annotate("text", 
               x = -Inf, 
               y = Inf, 
               label = gene_of_interest,
               hjust = -0.1, 
               vjust = 1.5, 
               size = 2,
               fontface = "bold",
               color = "white") + 
      theme(element_text(family = "Arial"))
  ) 
}, simplify = FALSE, USE.NAMES = TRUE)

combined_landmark_plot <- Reduce(`+`, landmark_plots) + 
  plot_layout(ncol=8) &
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = "black", color = NA),      # Entire plot background
        panel.background = element_rect(fill = "black", color = NA))

ggsave(filename = paste0("image/",
                         config$experiment_name,
                         "_combined_measured_landmarks.tiff"),
       combined_landmark_plot,
       device = "tiff",
       height = 11,
       width = 8.5,
       units = "in")

gene_of_interest <- "Mgp"

landmark_signal <- landmarks_count_matrix[gene_of_interest, ]

landmark_signal_table <- getSignalCoordinates(gene_of_interest,
                                              landmark_signal,
                                              height = height,
                                              width = width)

getSpatialPlot(
  landmark_signal_table,
  width = width,
  height = height,
  point_size = 0.5,
  channel_color_scale = setNames(c("red"), c(gene_of_interest)),
  background_file = paste0("image/",
                           config$experiment_name, 
                           "_background_scaled.jpg"),
  plot_legend = F
) 

# Plot imputed expression
cells_bins_probabilities <- Matrix::readMM(paste0("data/",
                                                  config$experiment_name,
                                                  "_cells_bins_probabilities.mtx"))

unique_bins_key <- read.csv(file = paste0("data/", 
                                          config$experiment_name,
                                          "_unique_bins_key.csv"))$index + 1

common_landmarks <- config$common_landmarks
# common_landmarks <- read.csv(file = paste0("data/", 
#                                            config$experiment_name,
#                                            "_kidney_landmarks_meta.csv"))$gene

gene_of_interest <- "Cspg4"
spatial_signal <- getSpatialSignal(as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
                                   cells_bins_probabilities,
                                   mapped.bins = NULL,
                                   log = F)

plotExpectedSignal(spatial_signal[unique_bins_key], 
                   gene_of_interest, 
                   height, 
                   width,
                   denoise = T)


p <- plotThresholdedSignal(
  spatial_signal[unique_bins_key],
  point_size = 0.5,
  gene_of_interest,
  height = height,
  width = width,
  background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
  percentile = 96,
  plot_legend = FALSE
) 
p

ggsave(filename = paste0("image/",
                         config$experiment_name,
                         "_",
                         stringr::str_to_lower(gene_of_interest),
                         ".tiff"),
       plot = p + theme(element_text(family = "Arial")),
       device = "tiff",
       height = 11,
       width = 8.5,
       units = "in")

common_nonzero_landmarks <- common_landmarks[rowSums(expression_matrix[common_landmarks, ]) != 0]

inferred_landmark_plots <- sapply(common_nonzero_landmarks, function(gene_of_interest) {
  spatial_signal <- getSpatialSignal(as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
                                     cells_bins_probabilities,
                                     mapped.bins = NULL,
                                     log = F)
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
               color = "white") + 
      theme(element_text(family = "Arial"))
  ) 
}, simplify = FALSE, USE.NAMES = TRUE)

combined_inferred_landmarks_plot <- Reduce(`+`, inferred_landmark_plots) + 
  plot_layout(ncol=8) &
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = "black", color = NA),      # Entire plot background
        panel.background = element_rect(fill = "black", color = NA))

ggsave(filename = paste0("image/",
                         config$experiment_name,
                         "_combined_inferred_landmarks.tiff"),
       combined_inferred_landmarks_plot,
       device = "tiff",
       height = 11,
       width = 8.5,
       units = "in")

interstitium_cluster_table <- read.csv(file = "../kidney_development/data/mouse_embryonic_interstitium_cluster_table.csv", 
                                       row.names = 1) %>% 
  dplyr::filter(stage == "e18.5d") %>%
  rownames_to_column(var = "cell_id") %>%
  mutate(cell_id = gsub("-e18.5d", "-wild type", cell_id)) %>%
  column_to_rownames(var = "cell_id")

cell_type_levels <- setNames(c("Endothelium",
                               "Interstitium",
                               "Leukocyte",
                               "Nephron epithelium",
                               "Podocyte",
                               "Ureteric epithelium",
                               "Stressed",
                               "Undetermined"),
                             levels(sce$cell_type_short))

sce$cell_type <- as.character(sce$cell_type)
sce$cell_type <- factor(cell_type_levels[as.character(sce$cell_type_short)], cell_type_levels)

cell_type_cluster_levels <- c(setdiff(levels(sce$cell_type), "Interstitium"), paste0("Interstitium_", c(paste0(1:7), paste0(8, "_", 1:2), "9", paste0(10, "_", 1:4), paste0(11:16))))
cell_type_clusters <- setNames(as.character(sce$cell_type), colnames(sce))
cell_type_clusters[row.names(interstitium_cluster_table)] <- paste0("Interstitium_", interstitium_cluster_table$cluster)
sce$cell_type_cluster <- factor(cell_type_clusters[colnames(sce)], levels = cell_type_cluster_levels)

cell_type_cluster <- "Interstitium_12"
spatial_signal <- getSpatialSignal(as.numeric(sce$cell_type_cluster == cell_type_cluster),
                                   cells_bins_probabilities,
                                   mapped.bins = NULL,
                                   log = F)



plotExpectedSignal(spatial_signal[unique_bins_key], 
                   cell_type_cluster, 
                   height, 
                   width,
                   denoise = T)


plotThresholdedSignal(
  spatial_signal[unique_bins_key],
  point_size = 0.5,
  cell_type_cluster,
  height = height,
  width = width,
  background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
  percentile = 99
)

cell_type_clusters <- paste0("Interstitium_", c(6, 14, "8_2", 12))

spatial_signals <- do.call(rbind, lapply(cell_type_clusters, function(cell_type_cluster) {
  return(getSpatialSignal(as.numeric(sce$cell_type_cluster == cell_type_cluster),
                          cells_bins_probabilities,
                          mapped.bins = NULL,
                          log = F))  
}))
row.names(spatial_signals) <- cell_type_clusters
spatial_signals <- spatial_signals[, unique_bins_key]

plotThresholdedSignals(spatial_signals,
                       height = config$scaled_height,
                       width = config$scaled_width,
                       percentiles = c(97, 97, 90, 95),
                       denoise = TRUE,
                       point_size = 1,
                       col = c("red", "white", "green", "yellow")[1:nrow(spatial_signals)],
                       plot_legend = FALSE)

ggsave(filename = paste0("image/mouse_",
                         stage, 
                         "_kidney_vascular_associated stroma.tiff"),
       device = "tiff",
       height = 8.5,
       width = 11,
       units = "in")
