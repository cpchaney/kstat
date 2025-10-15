#' Save spatial interaction plots (boxed and zoomed) to TIFF files.
#'
#' @param plots A named list of ggplot objects: `$boxed` and `$zoom`.
#' @param interaction_id Identifier string for the interaction being plotted.
#' @param config A list-like object containing configuration values (e.g., stage, image dimensions).
#' @param cell_types A character vector of the interacting cell types (used in file naming).
#' @param zoom_box A list or object with `x` and `y` vectors defining the zoomed region (pixel ranges).
#'
#' @return Saves plots to disk, no return value.

save_interaction_plots <- function(plots, interaction_id, config, cell_types, zoom_box) {
  # Construct output file prefix using config stage and cell type labels
  prefix <- file.path(
    "image",
    sprintf("mouse_%s_stat_%s_%s",
      config$stage,
      tolower(interaction_id),
      paste0(cell_types, collapse = "_"))
  )
  
  # Save the full (boxed) interaction plot
  ragg::agg_tiff(paste0(prefix, "_boxed.tiff"),
                 width = config$scaled_width,
                 height = config$scaled_height,
                 units = "px",
                 res = 300)
  print(plots$boxed)
  dev.off()
  
  # Compute dimensions for zoomed-in region
  zoom_width <- diff(zoom_box$x)
  zoom_height <- diff(zoom_box$y)

  # Save the zoomed-in plot (legend removed)
  ragg::agg_tiff(paste0(prefix, "_zoom.tiff"),
                 width = zoom_width,
                 height = zoom_height,
                 units = "px",
                 res = 300)
  print(plots$zoom + theme(legend.position = "none"))
  dev.off()
}
