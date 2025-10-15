library(tidyverse)
library(rjson)
library(here)
library(patchwork)
library(purrr)

# Source utility scripts
source(here("../tools/R/sc-rna-seq_util.R"))
source(here("R/utils.R"))
source(here("R/visualization.R"))

getSpatialPlot <- function(signal_table,
                           point_size,
                           channel_color_scale,
                           height,
                           width,
                           background_file = NULL,
                           plot_legend = TRUE, 
                           coord_flip = FALSE) {
  
  
  # Initialize ggplot object with background image and signal data table
  p <- ggplot(signal_table, aes(x = Y, y = X)) +
    coord_fixed(xlim = c(0, height), ylim = c(0, width), expand = FALSE)
  
  if (!is.null(background_file)) {
    # Load background image from file
    background_image <- jpeg::readJPEG(background_file)
    background_image <- as.cimg(background_image)
    background_image <- imrotate(background_image, angle = 90)
    background_image <- background_image[ , dim(background_image)[2]:1, , ]
    #background_image <- as.array(background_image)[, , 1, ]
    p <- p +
      ggpubr::background_image(background_image)
  }
  p <- p +
    geom_point(aes(colour = channel), 
               alpha = 0.8, 
               size = point_size, 
               stroke = 0) +
    scale_color_manual(values = channel_color_scale) +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'black'),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = -4, l = -4),
          legend.background = element_rect(fill = 'black'),
          legend.key = element_rect(fill = "black"),
          legend.text = element_text(size = 24, colour = 'white'),
          legend.title = element_blank(),
          legend.position = c(1, 0), 
          legend.justification = c(1, 0))
  
  # Hide the legend if plot_legend is set to FALSE
  if (!plot_legend) {
    p <- p + theme(legend.position = 'none')
  }
  
  return(p)
}

# Define age and load configuration
stage <- "p3d"
config <- here("config", str_c("mouse_", stage, ".json")) %>%
  fromJSON(file = .)

sce <- readRDS(file = paste0("data/mouse_",
                            stage,
                            "_snuc_sce.rds"))
expression_matrix <- assay(sce, "imputed")

width <- config$scaled_width
height <- config$scaled_height

# Load data
landmarks_count_matrix <- Matrix::readMM(paste0("data/mouse_", stage, "_landmarks_matrix.mtx"))
cells_bins_probabilities <- Matrix::readMM(paste0("data/",
                                                  config$experiment_name,
                                                  "_cells_bins_probabilities.mtx"))
unique_bins_key <- read.csv(file = paste0("data/", 
                                          config$experiment_name,
                                          "_unique_bins_key.csv"))$index + 1
common_landmarks <- config$common_landmarks

# Preprocess data
landmarks <- config$landmarks
row.names(landmarks_count_matrix) <- landmarks

#ordered_landmarks <- read.csv(file = "data/mouse_e18.5d_cartana_genes_ordered.csv")$gene
#ordered_landmarks <- c(ordered_landmarks, setdiff(landmarks, ordered_landmarks))
ordered_landmarks <- landmarks[order(landmarks)]

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
      plot_legend = FALSE,
      coord_flip = TRUE) +
      annotate("text", 
               x = height - 24,
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
  plot_layout(ncol=7) &
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

composite_width <- 7 * (height / 300)
composite_height <- 13 * (width / 300)
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
  units = "in",
  # Total height in inches
  dpi = 300,                        # Resolution in DPI
  bg = "transparent"
)

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
gene_of_interest <- "Smoc2"
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
  percentile = 31,
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
  threshold <- 100 * (1 - sum(landmarks_count_matrix[gene_of_interest, ]) / sum(spatial_signal[unique_bins_key] > 0))
  
  p <- plotThresholdedSignal(
    spatial_signal[unique_bins_key],
    gene_of_interest,
    height = height,
    width = width,
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    percentile = 97,
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
                         "_combined_og_stat_inferred_landmarks.tiff"),
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

# Figure generation
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

###############
## Scale bar ##
###############

rotate90_counterclockwise <- function(img) {
  if (length(dim(img)) == 3) {
    apply(img, 3, function(x) t(x)[ncol(x):1, ]) |> 
      array(dim = c(dim(img)[2], dim(img)[1], dim(img)[3]))
  } else {
    t(img)[ncol(img):1, ]
  }
}

rotate90_clockwise <- function(img) {
        if (length(dim(img)) == 3) {
                # Color image
                apply(img, 3, function(x) t(x)[, nrow(x):1]) |>
                        array(dim = c(dim(img)[2], dim(img)[1], dim(img)[3]))
        } else {
                # Grayscale image
                t(img)[, nrow(img):1]
        }
}

flip_horizontal <- function(img) {
  if (length(dim(img)) == 3) {
    # RGB image: flip each channel
    apply(img, 3, function(x) x[, ncol(x):1]) |>
      array(dim = dim(img))
  } else {
    # Grayscale image
    img[, ncol(img):1]
  }
}

background_image <- jpeg::readJPEG(paste0("image/", config$experiment_name, "_background_scaled.jpg"))
background_image <- flip_horizontal(background_image)
background_image <- rotate90_clockwise(background_image)
width <- ncol(background_image)
height <- nrow(background_image)

p <- ggplot() +
      ggpubr::background_image(background_image) +
  coord_fixed(
    xlim = c(0, width),
    ylim = c(0, height),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "black"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(
      t = 0,
      # Top margin
      r = 0,
      # Right margin
      b = 0,
      # Bottom margin
      l = 0,
      # Left margin
      unit = "mm"
    ),
    legend.position = "none"
  )

# === Setup: constants ===
scale_length_um <- 100
pixel_size_um <- 2.56
scale_length_units <- scale_length_um / pixel_size_um

bar_y <- 64 # vertical position of the bar (constant y)
bar_x_start <- 64 # where the bar starts
bar_x_end <- bar_x_start + scale_length_units
bar_height <- 2 # visual thickness in data units

p_scale_bar <- p +
  geom_rect(
    aes(
      xmin = bar_x_start,
      xmax = bar_x_end,
      ymin = bar_y - bar_height / 2,
      ymax = bar_y + bar_height / 2
    ),
    fill = "white"
  ) +
  annotate("text",
    x = (bar_x_start + bar_x_end) / 2,
    y = bar_y - 16, # adjust below bar
    label = paste0(scale_length_um, " µm"),
    size = 2,
    color = "white"
  )

ragg::agg_tiff(
  filename = paste0(
      "image/mouse_",
      config[["experiment_name"]],
      "_stat_rotated_scale_bar.tiff"
  ),
  width = width,
  height = height,
  units = "px",
  res = 300
)
p_scale_bar
dev.off()


################################
## snRNA-seq characterization ##
################################
coarse_cell_types <- as.character(sce$cell_type_short)
coarse_cell_types[coarse_cell_types == "Pod"] <- "NE"

sce$coarse_cell_type <- factor(coarse_cell_types)

if (length(grep("sling", colnames(colData(sce)))) > 0) {
  cell_table <- as(colData(sce)[, colnames(colData(sce))[-grep("sling", colnames(colData(sce)))]], "data.frame")
} else {
  cell_table <- as(colData(sce), "data.frame")
}
cell_table <- rownames_to_column(cell_table)
cell_table <- merge(x = cell_table, y = reducedDim(sce, "UMAP")[, 1:3], by.x = "rowname", by.y = "row.names")
colnames(cell_table)[(ncol(cell_table) - 2):ncol(cell_table)] <- paste0("UMAP_", 1:3)
cell_table <- merge(x = cell_table, y = reducedDim(sce, "DMAP"), by.x = "rowname", by.y = "row.names")
cell_table <- column_to_rownames(cell_table, "rowname")

cell_type_short_centers <- cell_table %>%
  group_by(cell_type_short) %>%
  dplyr::select(UMAP_1, UMAP_2, UMAP_3) %>%
  summarize_all(mean)

cell_type_centers <- cell_table %>%
  group_by(coarse_cell_type) %>%
  dplyr::select(UMAP_1, UMAP_2, UMAP_3) %>%
  summarize_all(mean)

expression_matrix <- assay(sce, "imputed")

cell_type_color_scale <- getFactorColors(factor(cell_table$coarse_cell_type))

ggplot(cell_table, aes(x = UMAP_1, y = UMAP_3)) +
        geom_point(aes(colour = coarse_cell_type), alpha = 1) +
        scale_color_manual(
                values = cell_type_color_scale,
                name = "Cell type"
        ) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 8))) +
        theme_bw() +
        theme(
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 16),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 20),
                axis.line = element_line(colour = "black")
        )

ggsave(
  filename = paste0("image/mouse_", stage, "_kidney_coarse_cell_type_umap_1_3.png"),
  width = 11,
  height = 8.5,
  units = "in"
)

##################################
## Heatmap of marker expression ##
##################################
## # EIGEN marker finding
## eigen_genes <- eigen(sce,
##   groups = "coarse_cell_type",
##   BPPARAM = BiocParallel::MulticoreParam(workers = parallel::detectCores() - 8)
## )

## dir_path <- "output" # specify your desired directory path
## # Create the directory if it doesn't exist
## if (!dir.exists(dir_path)) {
##   dir.create(dir_path, recursive = TRUE)
## }
## write.csv(eigen_genes,
##   file = "output/coarse_cell_type_eigen_genes.csv",
##   quote = FALSE
## )
## coarse_cell_type_markers <- eigen_genes[, -grep("p.value", colnames(eigen_genes))]

## genes_of_interest <- Reduce(union, lapply(levels(sce$coarse_cell_type), function(coarse_cell_type) {
##         row.names(coarse_cell_type_markers)[order(coarse_cell_type_markers[, coarse_cell_type])][1:6]
## }))

genes_of_interest <- intersect(
        read.csv("data/coarse_cell_type_markers.csv", header = FALSE)$V1,
        row.names(expression_matrix)
)

expr_df <- as.data.frame(as.matrix(t(expression_matrix[
  genes_of_interest,
  row.names(cell_table)
])))

meta_df <- rownames_to_column(cell_table, var = "cell_id") %>%
  dplyr::select("cell_id", "coarse_cell_type")

# 1. Add cell type to expression data
expr_long <- expr_df %>%
  rownames_to_column("cell_id") %>%
  pivot_longer(-cell_id, names_to = "gene", values_to = "expression") %>%
  inner_join(meta_df, by = "cell_id")

# 2. Compute mean expression per cell type
mean_expr <- expr_long %>%
  group_by(coarse_cell_type, gene) %>%
  summarize(mean_expr = mean(expression), .groups = "drop")

# 3. Pivot to wide format: genes x cell types
expr_matrix <- mean_expr %>%
  pivot_wider(names_from = coarse_cell_type, values_from = mean_expr) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# 4. Z-score scale each gene (row-wise)
scaled_expr <- t(scale(t(expr_matrix)))

# 5. Clip to range [-3, 3]
clipped_expr <- pmin(pmax(scaled_expr, -3), 3)

# 6. Plot heatmap
p <- pheatmap(clipped_expr,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = viridis(100),
  fontsize_row = 8,
  fontsize_col = 12,
  angle_col = 0,
  border_color = NA,
  main = "Marker Gene Expression by Cell Type (Z-scaled)"
)

png(paste0("image/_mouse_", stage, "_coarse_cell_type_marker_heatmap.png"), width = 2000, height = 1800, res = 300)
p
dev.off()



  
