library(tidyverse)
library(rjson)
library(here)
library(purrr)
library(rhdf5)
library(plotly)
library(ggplot2)
library(patchwork)
library(cowplot)

# Source utility scripts
source(here("../tools/R/sc-rna-seq_util.R"))
source(here("R/utils.R"))
source(here("R/visualization.R"))

use_condaenv("statlas", required = TRUE)

get_threshold <- function(input_vector, k = 1024) {
  # Filter for positive, unique values
  s <- sort(unique(input_vector[input_vector > 0]))
  n <- length(s)

  # Gracefully reduce k if too large
  max_k <- floor(n / 2) - 1
  if (k > max_k) {
    warning(sprintf("k = %d is too large for input (max valid k = %d). Using k = %d.", k, max_k, max_k))
    k <- max_k
  }

  if (n < 4) {
    warning("Not enough unique positive values to perform regression. Returning NA.")
    return(NA_real_)
  }

  m <- floor(n / 2)

  # Index ranges for two models
  x1 <- seq(m - k, m + k - 1)
  y1 <- s[x1]

  x2 <- seq(n - k + 1, n)
  y2 <- s[x2]

  # Fit two linear models
  model_1 <- lm(y1 ~ x1)
  model_2 <- lm(y2 ~ x2)

  # Extract coefficients
  a1 <- coef(model_1)["x1"]
  b1 <- coef(model_1)["(Intercept)"]

  a2 <- coef(model_2)["x2"]

  b2 <- coef(model_2)["(Intercept)"]

  # Compute intersection point
  if (a1 == a2) {
    warning("Lines are parallel; returning midpoint value.")
    x_val <- m
  } else {
    x_val <- (b2 - b1) / (a1 - a2)
  }

  # Round and clip index
  idx <- round(x_val)
  idx <- min(max(idx, 1), length(s))

  return(s[idx])
}

find_signal_knee <- function(signal, plot = FALSE) {
  # Step 1: Sort and get unique values
  signal <- sort(unique(signal))

  # Step 2: Construct x (index) and y (signal values)
  x <- seq_along(signal)
  y <- signal

  # Step 3: Line from first to last point
  point1 <- c(x[1], y[1])
  point2 <- c(x[length(x)], y[length(y)])

  # Step 4: Compute perpendicular distances to line
  line_vec <- point2 - point1
  line_vec_norm <- line_vec / sqrt(sum(line_vec^2))

  # Compute vector from point1 to each point
  vecs_to_line <- cbind(x - point1[1], y - point1[2])

  # Compute projection lengths onto line
  proj_lengths <- rowSums(vecs_to_line * matrix(rep(line_vec_norm, length(x)), ncol = 2, byrow = TRUE))

  # Projected points on the line
  proj_points <- matrix(rep(point1, length(x)), ncol = 2, byrow = TRUE) +
    proj_lengths %*% t(line_vec_norm)

  # Distances from actual points to their projection on the line
  distances <- sqrt(rowSums((cbind(x, y) - proj_points)^2))

  # Step 5: Find index of max distance ? the "knee"
  knee_index <- which.max(distances)
  threshold_value <- signal[knee_index]

  # Optional Plot
  if (plot) {
    plot(x, y, type = "l", main = "Knee Detection", xlab = "Index", ylab = "Signal")
    points(x[knee_index], y[knee_index], col = "red", pch = 19)
    lines(c(point1[1], point2[1]), c(point1[2], point2[2]), col = "blue", lty = 2)
    legend("bottomright",
      legend = c("Signal", "Knee", "Line"), col = c("black", "red", "blue"),
      lty = c(1, NA, 2), pch = c(NA, 19, NA)
    )
  }

  return(list(
    knee_index = knee_index,
    threshold_value = threshold_value
  ))
}


save_legend <- function(p, filename, dpi = 600, bg = "white") {
  stopifnot(
    requireNamespace("cowplot"),
    requireNamespace("grid"),
    requireNamespace("ragg")
  )

  # Get all guide-box components (handles multiple legends)
  legs <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
  if (length(legs) == 0) stop("No legend found. Is legend.position = 'none'?")

  # Pick the largest legend by area (in inches)
  areas <- vapply(legs, function(g) {
    w <- grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE)
    h <- grid::convertHeight(grid::grobHeight(g), "in", valueOnly = TRUE)
    w * h
  }, numeric(1))
  leg <- legs[[which.max(areas)]]

  w_in <- grid::convertWidth(grid::grobWidth(leg), "in", valueOnly = TRUE)
  h_in <- grid::convertHeight(grid::grobHeight(leg), "in", valueOnly = TRUE)

  # Choose device based on extension
  ext <- tolower(tools::file_ext(filename))
  device <- switch(ext,
    "tif"  = ragg::agg_tiff,
    "tiff" = ragg::agg_tiff,
    "png"  = ragg::agg_png,
    "svg"  = svglite::svglite,
    "pdf"  = grDevices::cairo_pdf,
    stop("Unsupported file extension: ", ext)
  )

  ggplot2::ggsave(
    filename,
    plot = cowplot::ggdraw(leg),
    width = w_in, height = h_in, units = "in",
    dpi = if (ext %in% c("png", "tif", "tiff")) dpi else NA,
    bg = bg,
    device = device
  )
}

# Define age and load configuration
stage <- "e18.5d"
config <- here("config", str_c("mouse_", stage, ".json")) %>%
  fromJSON(file = .)

sce <- readRDS(file = "../kidney_development/mouse_e18.5d/data/mouse_e18.5d_sn_sce.rds")
expression_matrix <- assay(sce, "imputed")
background_image <- jpeg::readJPEG(paste0("image/", config$experiment_name, "_background_scaled.jpg"))

width <- config$scaled_width
height <- config$scaled_height

# Plot imputed expression
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
# common_landmarks <- read.csv(file = paste0("data/",
#                                            config$experiment_name,
#                                            "_kidney_landmarks_meta.csv"))$gene

#################################
## Ligand receptor interaction ##
#################################
# interaction_id <- "CPI-SC045BF427E"
interaction_id <- "CPI-SC050D66C45"
# interaction_id <- "CPI-SC0675C8528"
# interaction_id <- "CPI-SC0B739B6BE"

interaction_df <- read.csv("data/ligand_receptor_components.csv")

cell_type_list <- paste0("NPC", 1:3)

hdf5_path <- paste0(
  "data/mouse_e18.5d_kidney_",
  paste0(cell_type_list, collapse = "_"),
  "_transport_distributions.h5"
)

# Construct dataset paths
ligand_path <- paste0("/", interaction_id, "/ligand")
delivered_path <- paste0("/", interaction_id, "/delivered_prop")
receptor_path <- paste0("/", interaction_id, "/receptor")
received_path <- paste0("/", interaction_id, "/received_prop")

# Read tensors
ligand <- as.numeric(h5read(hdf5_path, ligand_path))
ligand <- ligand / max(ligand)
delivered_prop <- as.numeric(h5read(hdf5_path, delivered_path))
receptor <- as.numeric(h5read(hdf5_path, receptor_path))
receptor <- receptor / max(receptor)
received_prop <- as.numeric(h5read(hdf5_path, received_path))
remaining_ligand <- as.numeric(delivered_prop != 0) - delivered_prop

ligand_channel_name <- dplyr::filter(interaction_df, id_cp_interaction == interaction_id & component_type == "ligand") %>%
  pull(mgi_symbol) %>%
  paste0(collapse = "+")

receptor_channel_name <- dplyr::filter(interaction_df, id_cp_interaction == interaction_id & component_type == "receptor") %>%
  pull(mgi_symbol) %>%
  paste0(collapse = "+")

ligand_mask <- numeric(length(ligand))
ligand_mask[which(ligand > quantile(ligand, 0.99))] <- 1

receptor_mask <- numeric(length(receptor))
receptor_mask[which(receptor > quantile(receptor, 0.99))] <- 1

ligand_signal_table <- getSignalCoordinates(ligand_channel_name, ligand * ligand_mask, height, width)

receptor_signal_table <- getSignalCoordinates(receptor_channel_name, receptor * receptor_mask, height, width)

signal_table <- rbind(
  ligand_signal_table,
  receptor_signal_table
)

channel_levels <- c(
  ligand_channel_name,
  receptor_channel_name,
  paste0(ligand_channel_name, "-", receptor_channel_name)
)
channel_color_scale <- setNames(
  c("red", "yellow", "white"),
  channel_levels
)

point_size <- 0.5

p <- ggplot(signal_table, aes(x = X, y = Y)) +
  ggpubr::background_image(background_image) +
  geom_point(aes(colour = channel, alpha = magnitude),
    size = point_size,
    stroke = 0
  ) +
  scale_color_manual(values = channel_color_scale) +
  scale_alpha(guide = "none", limits = c(0, 1)) +
  coord_fixed(
    xlim = c(0, dim(background_image)[2]),
    ylim = c(0, dim(background_image)[1]),
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

p <- p +
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

zoom_x <- c(582, 688)
zoom_y <- config$scaled_height - c(18, 112)

p_box <- p +
  expand_limits(x = zoom_x, y = zoom_y) +
  geom_rect(
    data = data.frame(
      xmin = zoom_x[1], xmax = zoom_x[2],
      ymin = zoom_y[2], ymax = zoom_y[1]
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "lightgrey", linewidth = 0.25, alpha = 1
  )

zoom_height <- zoom_y[1] - zoom_y[2]
zoom_width <- zoom_x[2] - zoom_x[1]

background_image_zoom <- background_image[zoom_x[1]:zoom_x[2], zoom_y[2]:zoom_y[1], ]

signal_table_zoom <- dplyr::filter(
  signal_table, X >= zoom_x[1] & X <= zoom_x[2]
) %>% dplyr::filter(Y >= zoom_y[2] & Y <= zoom_y[1])


p_zoom <- ggplot(signal_table_zoom, aes(x = X, y = Y)) +
  ggpubr::background_image(background_image_zoom) +
  geom_point(aes(colour = channel, alpha = magnitude),
    size = point_size,
    stroke = 0
  ) +
  scale_color_manual(values = channel_color_scale) +
  scale_alpha(guide = "none", limits = c(0, 1)) +
  coord_fixed(
    xlim = zoom_x,
    ylim = rev(zoom_y),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "black", colour = NA), # <- add this
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "black", colour = NA),
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
    legend.position = "none",
  )

ragg::agg_tiff(
  filename = paste0(
    "image/mouse_e18.5d_stat_",
    stringr::str_to_lower(
      paste0(
        interaction_id,
        "_",
        paste0(cell_type_list, collapse = "_")
      )
    ),
    "_ligand_receptor_boxed.tiff"
  ),
  width = config$scaled_width,
  height = config$scaled_height,
  units = "px",
  res = 300
)
p_box + theme(legend.position = "none")
dev.off()

p_zoom_bordered <- p_zoom +
  theme(
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 1.0),
    panel.background = element_rect(fill = "black", colour = NA)
  )

ragg::agg_tiff(
  filename = paste0(
    "image/mouse_e18.5d_stat_",
    stringr::str_to_lower(
      paste0(
        interaction_id,
        "_",
        paste0(cell_type_list, collapse = "_")
      )
    ),
    "_ligand_receptor_zoom.tiff"
  ),
  width = zoom_width,
  height = zoom_height,
  units = "px",
  res = 300
)
p_zoom_bordered + theme(legend.position = "none")
dev.off()

## save_legend(p + theme(legend.text = element_text(size = 8)),
##   filename = paste0(
##     "image/mouse_e18.5d_stat_",
##     stringr::str_to_lower(
##       paste0(
##         interaction_id,
##         "_",
##         paste0(cell_type_list, collapse = "_")
##       )
##     ),
##     "_transport_result_legend.tiff"
##   ),
##   dpi = 300
## )

#######################
## Transport results ##
#######################
delivered_ligand_signal_table <- getSignalCoordinates(
  ligand_channel_name,
  delivered_prop * ligand_mask,
  height,
  width
)

bound_receptor_signal_table <- getSignalCoordinates(
  paste0(ligand_channel_name, "-", receptor_channel_name),
  received_prop * receptor_mask,
  height,
  width
)

transport_signal_table <- rbind(
  delivered_ligand_signal_table,
  bound_receptor_signal_table
)

p <- ggplot(transport_signal_table, aes(x = X, y = Y)) +
  ggpubr::background_image(background_image) +
  geom_point(aes(colour = channel, alpha = magnitude),
    size = point_size,
    stroke = 0
  ) +
  scale_color_manual(values = channel_color_scale) +
  scale_alpha(guide = "none", limits = c(0, 1)) +
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

p <- p +
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

p_box <- p +
  expand_limits(x = zoom_x, y = zoom_y) +
  geom_rect(
    data = data.frame(
      xmin = zoom_x[1], xmax = zoom_x[2],
      ymin = zoom_y[2], ymax = zoom_y[1]
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, colour = "lightgrey", linewidth = 0.25, alpha = 1
  )

background_image_zoom <- background_image[zoom_x[1]:zoom_x[2], zoom_y[2]:zoom_y[1], ]

transport_signal_table_zoom <- dplyr::filter(
  transport_signal_table, X >= zoom_x[1] & X <= zoom_x[2]
) %>% dplyr::filter(Y >= zoom_y[2] & Y <= zoom_y[1])


p_zoom <- ggplot(transport_signal_table_zoom, aes(x = X, y = Y)) +
  ggpubr::background_image(background_image_zoom) +
  geom_point(aes(colour = channel, alpha = magnitude),
    size = point_size,
    stroke = 0
  ) +
  scale_color_manual(values = channel_color_scale) +
  scale_alpha(guide = "none", limits = c(0, 1)) +
  coord_fixed(
    xlim = zoom_x,
    ylim = rev(zoom_y),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "black", colour = NA), # <- add this
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "black", colour = NA),
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
    legend.position = "none",
  )

ragg::agg_tiff(
  filename = paste0(
    "image/mouse_e18.5d_stat_",
    stringr::str_to_lower(
      paste0(
        interaction_id,
        "_",
        paste0(cell_type_list, collapse = "_")
      )
    ),
    "_transport_result_boxed.tiff"
  ),
  width = config$scaled_width,
  height = config$scaled_height,
  units = "px",
  res = 300
)
p_box + theme(legend.position = "none")
dev.off()

p_zoom_bordered <- p_zoom +
  theme(
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 1.0),
    panel.background = element_rect(fill = "black", colour = NA)
  )

ragg::agg_tiff(
  filename = paste0(
    "image/mouse_e18.5d_stat_",
    stringr::str_to_lower(
      paste0(
        interaction_id,
        "_",
        paste0(cell_type_list, collapse = "_")
      )
    ),
    "_transport_result_zoom.tiff"
  ),
  width = zoom_width,
  height = zoom_height,
  units = "px",
  res = 300
)
p_zoom_bordered + theme(legend.position = "none")
dev.off()

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

p <- p +
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
    "image/mouse_e18.5d_stat_scale_bar.tiff"
  ),
  width = config$scaled_width,
  height = config$scaled_height,
  units = "px",
  res = 300
)
p
dev.off()


# ---- Load all function definitions ----
source(here("R/setup_environment.R"))
source(here("R/signal_threshold.R"))
source(here("R/load_interaction_data.R"))
source(here("R/compute_signals.R"))
source(here("R/plot_interaction.R"))
source(here("R/save_plots.R"))
source(here("R/save_legend.R")) # optional if you need legend export

# ---- Initialize environment ----
initialize_environment("statlas")

# ---- Load configuration ----
config <- fromJSON(file = here("config/mouse_e18.5d.json"))

# ---- Define inputs ----
cell_types <- paste0("NPC", 1:3)
interaction_id <- "CPI-SC0675C8528"
zoom_box <- list(x = c(732, 838), y = c(32, 128))

# ---- Generate and save plots ----
plots <- plot_interaction(interaction_id, config, cell_types, point_size = 0.2, zoom_box = zoom_box)
save_interaction_plots(plots, interaction_id, config, cell_types, zoom_box)
