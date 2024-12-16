library(ggplot2)
library(ggpubr)
library(jpeg)
library(imager)
library(viridis)
library(dplyr)

# Helper function for loading and transforming the background image
loadBackgroundImage <- function(background_file, coord_flip) {
  if (is.null(background_file)) return(NULL)

  background_image <- jpeg::readJPEG(background_file)

  if (coord_flip) {
    background_image <- as.cimg(background_image)
    background_image <- imrotate(background_image, angle = -90)
    background_image <- as.array(background_image)[, , 1, ]
  }

  return(background_image)
}

# Base function for creating a spatial plot
getSpatialPlot <- function(signal_table, point_size, channel_color_scale, height, width,
                           alpha = 0.5, background_file = NULL, plot_legend = TRUE,
                           coord_flip = FALSE) {
  # Initialize the base ggplot object
  aes_mapping <- if (coord_flip) aes(x = Y, y = width - X) else aes(x = X, y = Y)
  p <- ggplot(signal_table, aes_mapping) +
    coord_fixed(xlim = c(0, width), ylim = c(0, height), expand = FALSE)

  # Add background image if provided
  background_image <- loadBackgroundImage(background_file, coord_flip)
  if (!is.null(background_image)) {
    p <- p + ggpubr::background_image(background_image)
  }

  # Add signal points and customize plot appearance
  p <- p +
    geom_point(aes(colour = channel), alpha = alpha, size = point_size, stroke = 0) +
    scale_color_manual(values = channel_color_scale) +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = 'black'),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = -4, l = -4),
          legend.background = element_rect(fill = 'black'),
          legend.key = element_rect(fill = "black"),
          legend.text = element_text(size = 24, colour = 'white'),
          legend.title = element_blank(),
          legend.position = if (plot_legend) c(1, 0) else "none",
          legend.justification = c(1, 0))

  return(p)
}

# Function to threshold and plot a spatial signal
plotThresholdedSignal <- function(spatial_signal, signal_identifier, height, width,
                                  background_file = NULL, percentile = 95, denoise = TRUE,
                                  point_size = 0.5, alpha = 0.5, plot_legend = TRUE, col = 'red') {
  if (denoise) {
    spatial_signal <- gaussianBlur(spatial_signal, height = height, width = width)
  }

  # Apply thresholding
  nonzero_signal <- spatial_signal[spatial_signal > 0]
  threshold <- if (!is.null(percentile)) quantile(jitter(nonzero_signal), percentile / 100) else 0
  binary_signal <- as.integer(spatial_signal >= threshold)

  # Extract signal coordinates
  signal_table <- getSignalCoordinates(signal_identifier, binary_signal, height, width) %>%
    filter(magnitude > 0)

  channel_color_scale <- setNames(col, signal_identifier)

  return(getSpatialPlot(
    signal_table,
    point_size = point_size,
    channel_color_scale = channel_color_scale,
    height = height,
    width = width,
    alpha = alpha,
    background_file = background_file,
    plot_legend = plot_legend
  ))
}

# Function to plot multiple thresholded spatial signals
plotThresholdedSignals <- function(spatial_signals, height, width, percentiles = NULL,
                                   denoise = TRUE, alpha = 0.5, point_size = 0.5,
                                   col = c("red", "yellow", "white", "#E7298A"),
                                   plot_legend = TRUE) {
  if (is.null(percentiles)) percentiles <- rep(95, nrow(spatial_signals))
  signal_identifiers <- rownames(spatial_signals)
  names(percentiles) <- signal_identifiers

  # Optionally denoise signals
  if (denoise) {
    spatial_signals <- t(apply(spatial_signals, 1, function(signal) {
      gaussianBlur(signal, height, width)
    }))
  }

  # Calculate thresholds
  thresholds <- sapply(signal_identifiers, function(identifier) {
    spatial_signal <- spatial_signals[identifier, ]
    if (!is.null(percentiles)) {
      nonzero_signal <- spatial_signal[spatial_signal > 0]
      return(quantile(jitter(nonzero_signal), percentiles[identifier] / 100))
    } else {
      return(0)
    }
  })

  # Create binary signals
  binary_signals <- do.call(cbind, lapply(signal_identifiers, function(identifier) {
    as.integer(spatial_signals[identifier, ] >= thresholds[identifier])
  }))

  # Combine signal data into a single table
  signal_table <- do.call(rbind, lapply(seq_along(signal_identifiers), function(i) {
    getSignalCoordinates(signal_identifiers[i], binary_signals[, i], height, width) %>%
      filter(magnitude > 0)
  }))

  signal_table$channel <- factor(signal_table$channel, levels = signal_identifiers)
  channel_color_scale <- setNames(col, signal_identifiers)

  return(getSpatialPlot(
    signal_table,
    point_size = point_size,
    channel_color_scale = channel_color_scale,
    height = height,
    width = width,
    alpha = alpha,
    background_file = NULL,
    plot_legend = plot_legend
  ))
}

# Function to plot an expected signal as a raster
plotExpectedSignal <- function(spatial_signal, signal_identifier, height, width,
                               denoise = TRUE, blur_radius = 1, plot_legend = TRUE) {
  if (denoise) {
    spatial_signal <- gaussianBlur(spatial_signal, height = height, width = width, sigma = blur_radius)
  }

  # Extract signal coordinates
  signal_table <- getSignalCoordinates(signal_identifier, spatial_signal, height, width)

  # Create raster plot
  p <- ggplot(signal_table, aes(x = X, y = Y)) +
    geom_raster(aes(fill = magnitude), interpolate = TRUE) +
    scale_fill_viridis(name = signal_identifier, option = 'B') +
    coord_fixed(xlim = c(0, width), ylim = c(0, height), expand = FALSE) +
    theme_dark() +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'black'),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(t = 0, r = -1, b = -2, l = -2, unit = 'mm'),
      legend.position = if (plot_legend) c(1, 0) else "none",
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = 'black'),
      legend.title = element_text(colour = 'white', face = 'bold', size = 12),
      legend.text = element_text(colour = 'white', size = 10)
    )

  return(p)
}
