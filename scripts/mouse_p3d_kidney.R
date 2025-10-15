# ---- Setup & Config ----
library(tidyverse)
library(rjson)
library(here)
library(patchwork)
library(purrr)

# Load utility functions
source(here("../tools/R/sc-rna-seq_util.R"))
source(here("R/utils.R"))
source(here("R/visualization.R"))

# ---- Parameters ----
stage <- "p3d"
config <- fromJSON(file = here("config", paste0("mouse_", stage, ".json")))

# ---- Data Load & Preprocessing ----
sce <- readRDS(paste0("data/mouse_", stage, "_snuc_sce.rds"))
expression_matrix <- assay(sce, "imputed")

width <- config$scaled_width
height <- config$scaled_height

landmarks_count_matrix <- Matrix::readMM(paste0("data/mouse_", stage, "_landmarks_matrix.mtx"))
cells_bins_probabilities <- Matrix::readMM(paste0("data/", config$experiment_name, "_cells_bins_probabilities.mtx"))
unique_bins_key <- read.csv(paste0("data/", config$experiment_name, "_unique_bins_key.csv"))$index + 1
common_landmarks <- config$common_landmarks
landmarks <- config$landmarks
rownames(landmarks_count_matrix) <- landmarks
ordered_landmarks <- sort(landmarks)

# ---- Helper: Plotting Functions ----
plot_measured_landmark <- function(gene, width, height, point_size = 0.2) {
  signal <- landmarks_count_matrix[gene, ]
  coords <- getSignalCoordinates(gene, signal, height, width)

  getSpatialPlot(
    coords,
    width = width,
    height = height,
    point_size = point_size,
    channel_color_scale = setNames("red", gene),
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    plot_legend = FALSE,
    coord_flip = TRUE
  ) +
    annotate("text", x = height - 24, y = 24, label = gene, hjust = 1, vjust = 0, color = "white", size = 2, fontface = "bold")
}

plot_inferred_expression <- function(gene, percentile = 97, point_size = 0.2) {
  signal <- getSpatialSignal(as.numeric(expression_matrix[gene, , drop = FALSE]),
                             cells_bins_probabilities)

  plotThresholdedSignal(
    signal[unique_bins_key],
    signal_identifier = gene,
    height = height,
    width = width,
    background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"),
    percentile = percentile,
    point_size = point_size,
    plot_legend = FALSE
  ) +
    annotate("text", x = -Inf, y = Inf, label = gene, hjust = -0.1, vjust = 1.5, size = 2, fontface = "bold", color = "white") +
    theme(element_text(family = "Arial"))
}

save_plot <- function(plot, filename, width, height, scale_factor = 1, units = "in", dpi = 300) {
  ggsave(filename,
         plot = plot,
         width = width * scale_factor,
         height = height * scale_factor,
         dpi = dpi,
         units = units,
         device = "tiff",
         bg = "transparent")
}

# ---- Measured Landmark Composite Plot ----
landmark_plots <- map(ordered_landmarks, ~plot_measured_landmark(.x, width, height))

combined_plot <- wrap_plots(landmark_plots, ncol = 7) &
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black"),
    element_text(family = "Arial"),
    plot.margin = margin(0, 0, 0, 0)
  )

composite_width <- 7 * (height / 300)
composite_height <- 13 * (width / 300)
scale_factor <- 8.5 / composite_width

save_plot(
  combined_plot,
  paste0("image/", config$experiment_name, "_combined_measured_landmarks.tiff"),
  composite_width,
  composite_height,
  scale_factor
)

# ---- Individual Measured Gene Plot ----
plot_measured_landmark("Mgp", width, height, point_size = 0.5)

# ---- Inferred Expression: Single Gene ----
smoc2_plot <- plot_inferred_expression("Smoc2", percentile = 31, point_size = 0.5)
save_plot(smoc2_plot,
          paste0("image/", config$experiment_name, "_smoc2_inferred.tiff"),
          width = 8.5,
          height = 11)

# ---- Inferred Expression: Composite Landmark Panel ----
nonzero_landmarks <- common_landmarks[rowSums(expression_matrix[common_landmarks, ]) != 0]
inferred_plots <- map(nonzero_landmarks, ~plot_inferred_expression(.x))
combined_inferred <- wrap_plots(inferred_plots, ncol = 8) &
  theme(
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

save_plot(combined_inferred,
          paste0("image/", config$experiment_name, "_combined_inferred_landmarks.tiff"),
          width = 8.5,
          height = 11)

# ---- Cluster-Level Visualization ----
interstitium_clusters <- read.csv("../kidney_development/data/mouse_embryonic_interstitium_cluster_table.csv", row.names = 1) %>%
  filter(stage == "e18.5d") %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_id = gsub("-e18.5d", "-wild type", cell_id)) %>%
  column_to_rownames("cell_id")

# Add fine-grained clusters
cell_type_levels <- setNames(
  c("Endothelium", "Interstitium", "Leukocyte", "Nephron epithelium", "Podocyte", "Ureteric epithelium", "Stressed", "Undetermined"),
  levels(sce$cell_type_short)
)

sce$cell_type <- factor(cell_type_levels[as.character(sce$cell_type_short)], levels = cell_type_levels)
cell_type_clusters <- setNames(as.character(sce$cell_type), colnames(sce))
cell_type_clusters[rownames(interstitium_clusters)] <- paste0("Interstitium_", interstitium_clusters$cluster)
sce$cell_type_cluster <- factor(cell_type_clusters[colnames(sce)])

# Composite cluster-level plot
clusters_to_plot <- paste0("Interstitium_", c(6, 14, "8_2", 12))
spatial_signals <- map_dfr(clusters_to_plot, function(cluster) {
  sig <- getSpatialSignal(as.numeric(sce$cell_type_cluster == cluster), cells_bins_probabilities)
  sig <- sig[, unique_bins_key]
  rownames(sig) <- cluster
  sig
})

plotThresholdedSignals(spatial_signals,
                       height = config$scaled_height,
                       width = config$scaled_width,
                       percentiles = c(97, 97, 90, 95),
                       denoise = TRUE,
                       point_size = 1,
                       col = c("red", "white", "green", "yellow")[1:nrow(spatial_signals)],
                       plot_legend = FALSE)

ggsave(filename = paste0("image/mouse_", stage, "_kidney_vascular_associated_stroma.tiff"),
       device = "tiff",
       height = 8.5,
       width = 11,
       units = "in")

# ---- Supplemental Figure Panel ----
genes_of_interest <- c("Pou3f3", "Ass1", "Cldn11", "Apcdd1")
supplemental_plots <- map(genes_of_interest, ~plot_inferred_expression(.x, percentile = 99))
combined_supplemental <- wrap_plots(supplemental_plots, ncol = 4)

composite_width <- 4 * (height / 300)
composite_height <- 1 * (width / 300)
scale_factor <- 8.5 / composite_width

save_plot(
  combined_supplemental,
  paste0("image/", config$experiment_name, "_combined_inferred_", paste0(tolower(genes_of_interest), collapse = "_"), ".tiff"),
  composite_width,
  composite_height,
  scale_factor
)

# ---- Scale Bar ----
rotate90_clockwise <- function(img) {
  if (length(dim(img)) == 3) {
    apply(img, 3, function(x) t(x)[, nrow(x):1]) |>
      array(dim = c(dim(img)[2], dim(img)[1], dim(img)[3]))
  } else {
    t(img)[, nrow(img):1]
  }
}
flip_horizontal <- function(img) {
  if (length(dim(img)) == 3) {
    apply(img, 3, function(x) x[, ncol(x):1]) |>
      array(dim = dim(img))
  } else {
    img[, ncol(img):1]
  }
}

img <- jpeg::readJPEG(paste0("image/", config$experiment_name, "_background_scaled.jpg"))
img <- flip_horizontal(img) |> rotate90_clockwise()
width <- ncol(img); height <- nrow(img)

scale_length_um <- 100
pixel_size_um <- 2.56
scale_length_units <- scale_length_um / pixel_size_um

ggplot() +
  ggpubr::background_image(img) +
  coord_fixed(xlim = c(0, width), ylim = c(0, height), expand = FALSE) +
  geom_rect(aes(
    xmin = 64,
    xmax = 64 + scale_length_units,
    ymin = 64 - 1,
    ymax = 64 + 1
  ), fill = "white") +
  annotate("text", x = 64 + scale_length_units / 2, y = 48, label = "100 µm", color = "white", size = 2) +
  theme_void()

ragg::agg_tiff(
  filename = paste0("image/mouse_", config$experiment_name, "_stat_rotated_scale_bar.tiff"),
  width = width,
  height = height,
  units = "px",
  res = 300
)
dev.off()

# ---- snRNA-seq Visualization ----
coarse_cell_types <- as.character(sce$cell_type_short)
coarse_cell_types[coarse_cell_types == "Pod"] <- "NE"
sce$coarse_cell_type <- factor(coarse_cell_types)

cell_table <- rownames_to_column(as.data.frame(colData(sce)))
cell_table <- merge(cell_table, reducedDim(sce, "UMAP")[, 1:3], by.x = "rowname", by.y = "row.names")
colnames(cell_table)[(ncol(cell_table) - 2):ncol(cell_table)] <- paste0("UMAP_", 1:3)

color_scale <- getFactorColors(factor(cell_table$coarse_cell_type))

ggplot(cell_table, aes(x = UMAP_1, y = UMAP_3, color = coarse_cell_type)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = color_scale, name = "Cell type") +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_bw(base_size = 14)

ggsave(paste0("image/mouse_", stage, "_kidney_coarse_cell_type_umap_1_3.png"), width = 11, height = 8.5, units = "in")

# ---- Marker Heatmap ----
genes_of_interest <- intersect(
  read.csv("data/coarse_cell_type_markers.csv", header = FALSE)$V1,
  rownames(expression_matrix)
)

expr_df <- t(expression_matrix[genes_of_interest, ]) %>%
  as.data.frame() %>%
  rownames_to_column("cell_id")

meta_df <- rownames_to_column(cell_table, "cell_id") %>%
  select(cell_id, coarse_cell_type)

expr_long <- expr_df %>%
  pivot_longer(-cell_id, names_to = "gene", values_to = "expression") %>%
  inner_join(meta_df, by = "cell_id")

mean_expr <- expr_long %>%
  group_by(coarse_cell_type, gene) %>%
  summarize(mean_expr = mean(expression), .groups = "drop")

expr_matrix <- mean_expr %>%
  pivot_wider(names_from = coarse_cell_type, values_from = mean_expr) %>%
  column_to_rownames("gene") %>%
  as.matrix()

scaled_expr <- t(scale(t(expr_matrix)))
clipped_expr <- pmin(pmax(scaled_expr, -3), 3)

p <- pheatmap::pheatmap(clipped_expr,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = viridis::viridis(100),
  fontsize_row = 8,
  fontsize_col = 12,
  border_color = NA,
  main = "Marker Gene Expression by Cell Type (Z-scaled)"
)

png(paste0("image/_mouse_", stage, "_coarse_cell_type_marker_heatmap.png"), width = 2000, height = 180_
