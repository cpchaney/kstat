compute_signal_components <- function(data, width, height) {
  delivered_dist <- sparseMatrix(
    i = data$delivered_idx[2, ] + 1,
    j = data$delivered_idx[1, ] + 1,
    x = data$delivered,
    dims = c(height, width)
  )
  received_dist <- sparseMatrix(
    i = data$received_idx[2, ] + 1,
    j = data$received_idx[1, ] + 1,
    x = data$received,
    dims = c(height, width)
  )
  
  delivered <- as.vector(t(as.matrix(delivered_dist)))
  received <- as.vector(t(as.matrix(received_dist)))
  
  list(delivered = delivered, received = received)
}

plot_interaction <- function(interaction_id, config, cell_types, zoom_box = NULL) {
  # ---- Load tensors and derived data ----
  tensors <- load_interaction_data(interaction_id, config, cell_types)
  derived <- compute_signal_components(tensors, config$scaled_width, config$scaled_height)
  
  ligand <- tensors$ligand
  receptor <- tensors$receptor
  
  # ---- Thresholds ----
  ligand_threshold <- find_signal_knee(ligand)
  receptor_threshold <- find_signal_knee(receptor)
  
  # ---- Load background ----
  background_image <- jpeg::readJPEG(
    file.path("image", paste0(config$experiment_name, "_background_scaled.jpg"))
  )
  
  # ---- Extract and filter signals ----
  ligand_signal <- getSignalCoordinates("ligand", ligand, config$scaled_height, config$scaled_width) %>%
    filter(magnitude > ligand_threshold$threshold_value)
  
  receptor_signal <- getSignalCoordinates("receptor", receptor, config$scaled_height, config$scaled_width) %>%
    filter(magnitude > receptor_threshold$threshold_value)
  
  signal_table <- bind_rows(ligand_signal, receptor_signal)
  color_scale <- c("ligand" = "red", "receptor" = "yellow")
  
  # ---- Base plot ----
  p_base <- ggplot(signal_table, aes(X, Y)) +
    ggpubr::background_image(background_image) +
    geom_point(aes(colour = channel, alpha = magnitude), size = 0.2, stroke = 0) +
    scale_color_manual(values = color_scale) +
    scale_alpha(guide = "none", range = c(0.1, 1)) +
    coord_fixed(
      xlim = c(0, config$scaled_width),
      ylim = c(0, config$scaled_height),
      expand = FALSE
    ) +
    theme_void() +
    theme(plot.background = element_rect(fill = "black"))
  
  p_box <- p_base  # start as copy of base
  
  # ---- Add zoom rectangle if zoom_box is provided ----
  if (!is.null(zoom_box)) {
    zoom_x <- zoom_box$x
    zoom_y <- config$scaled_height - zoom_box$y
    
    p_box <- p_box +
      geom_rect(
        aes(xmin = zoom_x[1], xmax = zoom_x[2], ymin = zoom_y[2], ymax = zoom_y[1]),
        fill = NA, colour = "lightgrey", linewidth = 0.25
      )
    
    # ---- Zoom plot ----
    zoomed_image <- background_image[zoom_x[1]:zoom_x[2], zoom_y[2]:zoom_y[1], , drop = FALSE]
    signal_zoom <- signal_table %>%
      filter(between(X, zoom_x[1], zoom_x[2]), between(Y, zoom_y[2], zoom_y[1]))
    
    p_zoom <- ggplot(signal_zoom, aes(X, Y)) +
      ggpubr::background_image(zoomed_image) +
      geom_point(aes(colour = channel, alpha = magnitude), size = 0.2) +
      scale_color_manual(values = color_scale) +
      coord_fixed(xlim = zoom_x, ylim = rev(zoom_y), expand = FALSE) +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "black"),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.25)
      )
  } else {
    p_zoom <- NULL
  }
  
  # ---- Return plots ----
  return(list(full = p_base, boxed = p_box, zoom = p_zoom))
}
