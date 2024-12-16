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

# Load the data
landmarks_count_matrix <- Matrix::readMM(paste0("data/mouse_", stage, "_landmarks_matrix.mtx"))
cells_bins_probabilities <- Matrix::readMM(paste0("data/",
                                                  config$experiment_name,
                                                  "_cells_bins_probabilities.mtx"))
unique_bins_key <- read.csv(file = paste0("data/", 
                                          config$experiment_name,
                                          "_unique_bins_key.csv"))$index + 1
common_landmarks <- config$common_landmarks

# Prepare data
landmarks <- config$landmarks
row.names(landmarks_count_matrix) <- landmarks

ordered_landmarks <- read.csv(file = "data/mouse_e18.5d_cartana_genes_ordered.csv")$gene
ordered_landmarks <- c(ordered_landmarks, setdiff(landmarks, ordered_landmarks))

# Plot measured expression
landmark_plots <- sapply(ordered_landmarks, function(gene_of_interest) {
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
      plot_legend = FALSE) +
      annotate("text", 
               x = width - 24,
               y = 0 + 24, 
               label = gene_of_interest, 
               hjust = 1, 
               vjust = 0,
               color = "white",
               size = 2,
               fontface = "bold")
    )
}, simplify = FALSE, USE.NAMES = TRUE)

combined_landmark_plot <- Reduce(`+`, landmark_plots) +
  plot_layout(ncol=15) &
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
    )) 

#combined_landmark_plot <- patchwork::wrap_plots(landmark_plots, ncol = 15)

# tiff(filename = paste0("image/",
#                          config$experiment_name,
#                          "_combined_measured_landmarks.tiff"),
#        height = 1080,
#        width = 1920,
#        type = "cairo")
# combined_landmark_plot
# dev.off()

composite_width <- 15 * (width / 300)
composite_height <- 6 * (height / 300)
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
  dpi = 300,                        # Resolution in DPI
  bg = "transparent"
)

gene_of_interest <- "H19"

landmark_signal <- landmarks_count_matrix[gene_of_interest, ]

landmark_signal_table <- getSignalCoordinates(gene_of_interest,
                                              landmark_signal,
                                              height = height,
                                              width = width)

p <- getSpatialPlot(
  landmark_signal_table,
  width = width,
  height = height,
  point_size = 0.1,
  channel_color_scale = setNames(c("red"), c(gene_of_interest)),
  background_file = paste0("image/",
                           config$experiment_name, 
                           "_background_scaled.jpg"),
  plot_legend = F
) 
p
ggsave(filename = paste0("image/",
                         config$experiment_name,
                         "_",
                         stringr::str_to_lower(gene_of_interest),
                         "_measured.tiff"),
       plot = p,
       width = width,
       height = height,
       units = "px",
       bg = "transparent",
       dpi = 300
)

# Plot imputed expression
gene_of_interest <- "Gucy1a1"
spatial_signal <- getSpatialSignal(as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
                                   cells_bins_probabilities,
                                   mapped.bins = NULL,
                                   log = F)

p <- plotExpectedSignal(spatial_signal[unique_bins_key], 
                   gene_of_interest, 
                   height, 
                   width,
                   denoise = T,
                   plot_legend = FALSE)
p
ggsave(filename = paste0("image/",
                       config$experiment_name,
                       "_",
                       stringr::str_to_lower(gene_of_interest),
                       "_expected.tiff"),
       plot = p,
       width = width,
       height = height,
       units = "px",
       bg = "transparent",
       dpi = 300
)

threshold <- 99
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
ggsave(filename = paste0("image/",
                       config$experiment_name,
                       "_",
                       stringr::str_to_lower(gene_of_interest),
                       "_inferred_",
                       threshold,
                       ".tiff"),
       plot = p,
       width = width,
       height = height,
       units = "px",
       bg = "transparent",
       dpi = 300
)

# Plot inferred expression of multiple genes simultaneously
genes_of_interest <- c("Gucy1a2", "Clca3a1", "Cldn11")
spatial_signals <- do.call(rbind, lapply(genes_of_interest, function(gene_of_interest) {
  return(getSpatialSignal(expression_matrix[gene_of_interest, ],
                          cells_bins_probabilities,
                          mapped.bins = NULL,
                          log = F))  
}))
row.names(spatial_signals) <- genes_of_interest
spatial_signals <- spatial_signals[, unique_bins_key]

plotThresholdedSignals(spatial_signals,
                       height = config$scaled_height,
                       width = config$scaled_width,
                       percentiles = c(98, 98, 98, 98),
                       denoise = TRUE,
                       alpha = 1,
                       point_size = 0.1025262,
                       col = c("white", "red", "yellow")[1:nrow(spatial_signals)],
                       plot_legend = FALSE)

ggsave(filename = paste0("image/mouse_",
                         stage, 
                         "_kidney_combined_inferred_",
                         paste0(stringr::str_to_lower(genes_of_interest), collapse = "_"),
                         ".tiff"),
       device = "tiff",
       height = height,
       width = width,
       units = "px",
       dpi = 300)


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

# Pericyte analysis
# interstitium_cluster_table <- read.csv(file = "../kidney_development/data/mouse_embryonic_interstitium_cluster_table.csv", 
#                                        row.names = 1) %>% 
#   dplyr::filter(stage == "e18.5d") %>%
#   rownames_to_column(var = "cell_id") %>%
#   mutate(cell_id = gsub("-e18.5d", "-wild type", cell_id)) %>%
#   column_to_rownames(var = "cell_id")

#interstitium_cluster_table[interstitium_cluster_table$cluster %in% c("8_1", "8_2"), "cluster"] <- "8"

# cell_type_levels <- setNames(c("Endothelium",
#                                "Interstitium",
#                                "Leukocyte",
#                                "Nephron epithelium",
#                                "Podocyte",
#                                "Ureteric epithelium",
#                                "Stressed",
#                                "Undetermined"),
#                              levels(sce$cell_type_short))


cell_type_cluster <- "Int_8_2"
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
  point_size = 1,
  cell_type_cluster,
  height = height,
  width = width,
  col = c("green"),
  background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
  percentile = 75,
  plot_legend = FALSE
)

ggsave(filename = paste0("image/mouse_",
                         stage, 
                         "_kidney_", 
                         stringr::str_to_lower(cell_type_cluster), 
                         ".tiff"),
       device = "tiff",
       height = height,
       width = width,
       units = "px",
       dpi = 300)

# Pericytes
cell_type_clusters <- paste0("Int_", c(6, 14, "8_2", 12))
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
                       percentiles = c(10, 5, 75, 85),
                       denoise = TRUE,
                       point_size = 0.5,
                       col = c("red", "white", "green", "yellow")[1:nrow(spatial_signals)],
                       plot_legend = FALSE)

ggsave(filename = paste0("image/mouse_",
                         stage, 
                         "_kidney_vascular_associated_stroma.tiff"),
       device = "tiff",
       height = height,
       width = width,
       units = "px",
       dpi = 300)

# Proximal tubule stroma
cell_type_clusters <- paste0("Int_", c(6, 7, 8))

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
                       percentiles = c(9, 95, 95),
                       denoise = TRUE,
                       point_size = 0.5,
                       col = c("red", "white", "yellow")[1:nrow(spatial_signals)],
                       plot_legend = FALSE)

ggsave(filename = paste0("image/mouse_",
                         stage, 
                         "_kidney_proximal_tubule_stroma.tiff"),
       device = "tiff",
       height = height,
       width = width,
       units = "px",
       dpi = 300)

# Proximal tubule and proximal tubule stroma
proximal_tubule <- as.numeric(sce$cell_type_cluster %in% paste0("PT", 1:3))
proximal_tubule_stroma <- as.numeric(sce$cell_type_cluster %in% paste0("Int_", c(6:7, paste0("8_", 1:2))))

spatial_signals <- rbind(getSpatialSignal(proximal_tubule,
                                          cells_bins_probabilities,
                                          mapped.bins = NULL,
                                          log = F),
                         getSpatialSignal(proximal_tubule_stroma,
                                          cells_bins_probabilities,
                                          mapped.bins = NULL,
                                          log = F))

row.names(spatial_signals) <- c("PT", "PTS")
spatial_signals <- spatial_signals[, unique_bins_key]

plotThresholdedSignals(spatial_signals,
                       height = config$scaled_height,
                       width = config$scaled_width,
                       percentiles = c(90, 95),
                       denoise = TRUE,
                       point_size = 0.5,
                       col = c("red", "yellow")[1:nrow(spatial_signals)],
                       plot_legend = FALSE)

ggsave(filename = paste0("image/mouse_",
                         stage, 
                         "_kidney_proximal_tubule_proximal_tubule_stroma.tiff"),
       device = "tiff",
       height = height,
       width = width,
       units = "px",
       dpi = 300)

# Figure 2 generation
genes_of_interest <- c("Podxl", "Smoc2", "Slc27a2", "Egfl7")
gene_panels <- sapply(genes_of_interest, function(gene_of_interest) {
  panel_plots <- list()
  landmark_signal <- landmarks_count_matrix[gene_of_interest, ]
  
  landmark_signal_table <- getSignalCoordinates(gene_of_interest,
                                                landmark_signal,
                                                height = height,
                                                width = width)
  
  p <- getSpatialPlot(
    landmark_signal_table,
    width = width,
    height = height,
    point_size = 0.1,
    channel_color_scale = setNames(c("red"), c(gene_of_interest)),
    background_file = paste0("image/",
                             config$experiment_name, 
                             "_background_scaled.jpg"),
    plot_legend = F
  ) 
  panel_plots <- c(panel_plots, list(p))
  
  spatial_signal <- getSpatialSignal(as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]),
                                     cells_bins_probabilities,
                                     mapped.bins = NULL,
                                     log = F)
  
  p <- plotExpectedSignal(spatial_signal[unique_bins_key], 
                          gene_of_interest, 
                          height, 
                          width,
                          denoise = T,
                          plot_legend = FALSE)
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
    plot_layout(ncol=1, widths=1) &
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.background = element_rect(fill = "black", color = NA),      # Entire plot background
          panel.background = element_rect(fill = "black", color = NA))
  
  ggsave(filename = paste0("image/",
                           config$experiment_name,
                           "_",
                           stringr::str_to_lower(gene_of_interest),
                           "_paneled.tiff"),
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
  ) +
    annotate("text", 
             x = height - 24,
             y = 0 + 24, 
             label = gene_of_interest, 
             hjust = 1, 
             vjust = 0,
             color = "white",
             size = 4,
             fontface = "bold")
  
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
    )) 

composite_width <- 4 * (height / 300)
composite_height <- 1 * (width / 300)
scale_factor <- 8.5 / composite_width

ggsave(
  filename = paste0("image/",
                    config$experiment_name,
                    "_combined_inferred_",
                    paste0(stringr::str_to_lower(genes_of_interest), collapse = "_"),
                    ".tiff"),
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
  dpi = 300,                        # Resolution in DPI
  bg = "transparent"
)

cluster_means<- cbind(as.data.frame(as.matrix(t(expression_matrix[common_landmarks, ]))), 
                           as.data.frame(colData(sce)[, "cluster", drop = FALSE])) %>%
  group_by(cluster) %>%
  summarize_all(mean) %>%
  column_to_rownames(var = "cluster") 

landmark_variances <- setNames(colVars(as.matrix(cluster_means)), common_landmarks)
landmark_variances[order(landmark_variances)]

landmark_means <- setNames(colMeans(as.matrix(cluster_means)), common_landmarks)
landmark_means[order(landmark_means, decreasing = TRUE)]

landmark_lower_quartiles <- setNames(rowQuantiles(expression_matrix[common_landmarks, ], probs=0.25), common_landmarks)
landmark_lower_quartiles[order(landmark_lower_quartiles, decreasing = TRUE)]

landmark_table <- data.frame(landmark = common_landmarks,
                             mean = landmark_means,
                             cluster_variance = landmark_variances)
                            
ggplot(landmark_table, aes(x = mean, y = cluster_variance)) +
  geom_point(aes(color = ifelse(landmark %in% c("H19", "Fxyd2", "Acta2", "Ldhb", "Cdkn1c"), "red", "black")), size = 2) +
  geom_text(
    aes(label = ifelse(landmark %in% c("H19", "Fxyd2", "Acta2", "Ldhb", "Cdkn1c"), landmark, "")),  # Add labels for points meeting the condition
    vjust = -1                                    # Adjust vertical position of labels
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
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"))

  
