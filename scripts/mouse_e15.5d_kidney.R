# ---- Setup & Configuration ----
library(tidyverse)
library(rjson)
library(here)
library(patchwork)
library(purrr)
library(pheatmap)

# Load custom utilities
source(here("../tools/R/sc-rna-seq_util.R"))
source(here("R/utils.R"))
source(here("R/visualization.R"))

# ---- Helper: Thresholded Spatial Plot ----
plotThresholdedSignal <- function(spatial_signal,
                                  signal_identifier,
                                  height,
                                  width,
                                  background_file = NULL,
                                  percentile = 95,
                                  denoise = TRUE,
                                  point_size = 0.5,
                                  plot_legend = TRUE,
                                  col = "red",
                                  coord_flip = FALSE) {
  if (denoise) {
    spatial_signal <- gaussianBlur(spatial_signal, height = height, width = width)
  }
  threshold <- if (!is.null(percentile)) {
    nonzero <- jitter(spatial_signal[spatial_signal > 0])
    quantile(nonzero, percentile / 100)
  } else 0

  binary_signal <- as.integer(spatial_signal >= threshold)
  signal_table <- getSignalCoordinates(signal_identifier, binary_signal, height, width)
  signal_table <- signal_table[signal_table$magnitude > 0, ]
  color_scale <- setNames(col, signal_identifier)

  getSpatialPlot(signal_table,
    point_size = point_size,
    channel_color_scale = color_scale,
    height = height,
    width = width,
    background_file = background_file,
    plot_legend = plot_legend,
    coord_flip = coord_flip
  )
}

# ---- Data Loading & Preprocessing ----
stage <- "e15.5d"
config <- fromJSON(file = here("config", str_c("mouse_", stage, ".json")))

sce <- readRDS(file = paste0("../kidney_development/data/mouse_", stage, "_sce.rds"))
expression_matrix <- assay(sce, "imputed")

width <- config$scaled_width
height <- config$scaled_height

# Load matrices
landmarks_count_matrix <- Matrix::readMM(paste0("data/mouse_", stage, "_landmarks_matrix.mtx"))
cells_bins_probabilities <- Matrix::readMM(paste0("data/", config$experiment_name, "_cells_bins_probabilities.mtx"))
unique_bins_key <- read.csv(paste0("data/", config$experiment_name, "_unique_bins_key.csv"))$index + 1
common_landmarks <- config$common_landmarks
landmarks <- config$landmarks
row.names(landmarks_count_matrix) <- landmarks
ordered_landmarks <- sort(landmarks)

# ---- Utility Plot Functions ----
plot_measured_landmark <- function(gene, width, height, flip = TRUE) {
  signal <- landmarks_count_matrix[gene, ]
  table <- getSignalCoordinates(gene, signal, height, width)
  getSpatialPlot(
    table, width, height,
    point_size = 0.5,
    channel_color_scale = setNames("red", gene),
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    plot_legend = FALSE,
    coord_flip = flip
  ) +
    annotate("text",
      x = width - 24, y = 24,
      label = gene, color = "white", hjust = 1, vjust = 0,
      size = 2, fontface = "bold"
    )
}

save_plot <- function(plot, filename, width, height, scale_factor = 1, dpi = 300, units = "in") {
  ggsave(filename,
    plot = plot,
    device = "tiff",
    width = width * scale_factor,
    height = height * scale_factor,
    dpi = dpi,
    units = units,
    bg = "transparent"
  )
}

# ---- Measured Landmark Composite ----
landmark_plots <- map(ordered_landmarks, ~plot_measured_landmark(.x, width, height))
n_col <- 10

combined_landmark_plot <- wrap_plots(landmark_plots, ncol = n_col) &
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black"),
    element_text(family = "Arial"),
    plot.margin = margin(0, 0, 0, 0)
  )

composite_width <- n_col * (width / 300)
composite_height <- ceiling(length(landmarks) / n_col) * (height / 300)
scale_factor <- 8.5 / composite_width

save_plot(
  combined_landmark_plot,
  paste0("image/", config$experiment_name, "_combined_measured_landmarks.tiff"),
  composite_width, composite_height,
  scale_factor
)

# ---- Individual Measured Landmark ----
plot_measured_landmark("Mgp", width, height)

# ---- Inferred (Imputed) Expression Example ----
plot_inferred_expression <- function(gene, percentile = 95) {
  signal <- getSpatialSignal(as.numeric(expression_matrix[gene, , drop = FALSE]),
    cells_bins_probabilities, mapped.bins = NULL, log = FALSE
  )
  plotThresholdedSignal(signal[unique_bins_key],
    signal_identifier = gene,
    height = height,
    width = width,
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    percentile = percentile,
    plot_legend = FALSE
  ) +
    annotate("text",
      x = height - 24, y = 24, label = gene,
      color = "white", size = 4, fontface = "bold"
    )
}

# Example gene: Smoc2
smoc2_plot <- plot_inferred_expression("Smoc2", percentile = 31)
save_plot(smoc2_plot,
  paste0("image/", config$experiment_name, "_smoc2_inferred.tiff"),
  width = 8.5, height = 11
)

# ---- Combined Inferred Landmark Plot ----
common_nonzero <- common_landmarks[rowSums(expression_matrix[common_landmarks, ]) != 0]
inferred_plots <- map(common_nonzero, ~plot_inferred_expression(.x, percentile = 97))
combined_inferred <- wrap_plots(inferred_plots, ncol = 8) &
  theme(
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

save_plot(
  combined_inferred,
  paste0("image/", config$experiment_name, "_combined_inferred_landmarks.tiff"),
  width = 8.5, height = 11
)

# ---- Interstitium Cluster Visualization ----
interstitium_clusters <- read.csv("../kidney_development/data/mouse_embryonic_interstitium_cluster_table.csv", row.names = 1) %>%
  filter(stage == "e18.5d") %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_id = gsub("-e18.5d", "-wild type", cell_id)) %>%
  column_to_rownames("cell_id")

cell_type_levels <- setNames(
  c("Endothelium", "Interstitium", "Leukocyte", "Nephron epithelium", "Podocyte",
    "Ureteric epithelium", "Stressed", "Undetermined"),
  levels(sce$cell_type_short)
)

sce$cell_type <- factor(cell_type_levels[as.character(sce$cell_type_short)], levels = cell_type_levels)
cell_type_clusters <- setNames(as.character(sce$cell_type), colnames(sce))
cell_type_clusters[row.names(interstitium_clusters)] <- paste0("Interstitium_", interstitium_clusters$cluster)
sce$cell_type_cluster <- factor(cell_type_clusters[colnames(sce)])

# Example cluster visualization
cluster <- "Interstitium_12"
signal <- getSpatialSignal(as.numeric(sce$cell_type_cluster == cluster), cells_bins_probabilities)
plotThresholdedSignal(signal[unique_bins_key], cluster, height, width,
  paste0("image/", config$experiment_name, "_background_scaled.jpg"), percentile = 99
)

# ---- Supplemental Figure Generation ----
genes_of_interest <- c("Pou3f3", "Ass1", "Cldn11", "Apcdd1")
supplemental_plots <- map(genes_of_interest, ~plot_inferred_expression(.x, percentile = 99))
combined_supplemental <- wrap_plots(supplemental_plots, ncol = 4)

save_plot(
  combined_supplemental,
  paste0("image/", config$experiment_name, "_combined_inferred_", paste0(str_to_lower(genes_of_interest), collapse = "_"), ".tiff"),
  width = 4 * (height / 300),
  height = 1 * (width / 300),
  scale_factor = 8.5 / (4 * (height / 300))
)

# ---- Scale Bar Generation ----
rotate90_clockwise <- function(img) {
  if (length(dim(img)) == 3) {
    apply(img, 3, function(x) t(x)[, nrow(x):1]) |> array(dim = c(dim(img)[2], dim(img)[1], dim(img)[3]))
  } else t(img)[, nrow(img):1]
}

background <- jpeg::readJPEG(paste0("image/", config$experiment_name, "_background_scaled.jpg"))
background <- rotate90_clockwise(background)
width <- ncol(background); height <- nrow(background)

scale_bar <- ggplot() +
  ggpubr::background_image(background) +
  coord_fixed(xlim = c(0, width), ylim = c(0, height), expand = FALSE) +
  theme_void() +
  geom_rect(aes(xmin = 48, xmax = 48 + 100 / 2.56, ymin = 47, ymax = 49), fill = "white") +
  annotate("text", x = 98, y = 32, label = "100 µm", color = "white", size = 2)

ragg::agg_tiff(
  filename = paste0("image/mouse_", config$experiment_name, "_stat_scale_bar.tiff"),
  width = width, height = height, units = "px", res = 300
)
scale_bar
dev.off()

# ---- snRNA-seq Characterization ----
cell_table <- as.data.frame(colData(sce))
cell_table <- rownames_to_column(cell_table)
cell_table <- merge(cell_table, reducedDim(sce, "UMAP")[, 1:3], by.x = "rowname", by.y = "row.names")
colnames(cell_table)[(ncol(cell_table) - 2):ncol(cell_table)] <- paste0("UMAP_", 1:3)

color_scale <- getFactorColors(factor(cell_table$coarse_cell_type))

ggplot(cell_table, aes(x = UMAP_1, y = UMAP_3, color = coarse_cell_type)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = color_scale, name = "Cell type") +
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(), legend.position = "right")

ggsave(paste0("image/mouse_", stage, "_kidney_coarse_cell_type_umap_1_3.png"), width = 11, height = 8.5, units = "in")

# ---- Marker Gene Heatmap ----
genes_of_interest <- intersect(
  read.csv("data/coarse_cell_type_markers.csv", header = FALSE)$V1,
  row.names(expression_matrix)
)

expr_long <- expression_matrix[genes_of_interest, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  pivot_longer(-cell_id, names_to = "gene", values_to = "expression") %>%
  inner_join(cell_table %>% select(rowname, coarse_cell_type), by = c("cell_id" = "rowname"))

expr_mean <- expr_long %>%
  group_by(coarse_cell_type, gene) %>%
  summarize(mean_expr = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = coarse_cell_type, values_from = mean_expr) %>%
  column_to_rownames("gene") %>%
  as.matrix()

scaled_expr <- t(scale(t(expr_mean)))
scaled_expr <- pmin(pmax(scaled_expr, -3), 3)

p <- pheatmap(scaled_expr,
  cluster_rows = TRUE, cluster_cols = FALSE,
  color = viridis::viridis(100), border_color = NA,
  main = "Marker Gene Expression by Cell Type (Z-scaled)"
)

png(paste0("image/mouse_", stage, "_coarse_cell_type_marker_heatmap.png"), width = 2000, height = 1800, res = 300)
p
dev.off()
