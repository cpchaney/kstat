#' Save the legend of a ggplot as a standalone image
#'
#' @param p A ggplot object containing a legend.
#' @param filename File name (including extension: .png, .tiff, .pdf, .svg).
#' @param dpi Dots per inch for raster formats (default = 600).
#' @param bg Background color (default = "white").
#'
#' @return Invisibly returns the legend grob.
#' @export
#'
#' @examples
#' p <- ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) + geom_point()
#' save_legend(p, "legend.png")
save_legend <- function(p, filename, dpi = 600, bg = "white") {
  stopifnot(inherits(p, "ggplot"))

  # Load required packages quietly
  required_pkgs <- c("cowplot", "grid", "ragg", "ggplot2")
  missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "))
  }

  # Extract all guide-box components (handles multiple legends)
  legs <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
  if (length(legs) == 0) stop("No legend found in the provided ggplot (legend.position may be 'none').")

  # Choose the largest legend by area (useful if multiple)
  areas <- vapply(legs, function(g) {
    w <- grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE)
    h <- grid::convertHeight(grid::grobHeight(g), "in", valueOnly = TRUE)
    w * h
  }, numeric(1))
  leg <- legs[[which.max(areas)]]

  # Get legend size in inches
  w_in <- grid::convertWidth(grid::grobWidth(leg), "in", valueOnly = TRUE)
  h_in <- grid::convertHeight(grid::grobHeight(leg), "in", valueOnly = TRUE)

  # Choose graphics device by file extension
  ext <- tolower(tools::file_ext(filename))
  device <- switch(ext,
    "png"  = ragg::agg_png,
    "tif"  = ragg::agg_tiff,
    "tiff" = ragg::agg_tiff,
    "svg"  = svglite::svglite,
    "pdf"  = grDevices::cairo_pdf,
    stop("Unsupported file extension: ", ext)
  )

  # Save the legend using ggsave
  ggplot2::ggsave(
    filename = filename,
    plot = cowplot::ggdraw(leg),
    width = w_in,
    height = h_in,
    units = "in",
    dpi = if (ext %in% c("png", "tif", "tiff")) dpi else NA,
    bg = bg,
    device = device
  )

  invisible(leg)
}
