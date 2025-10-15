#' Detect the 'knee' (elbow) point in a sorted signal curve.
#'
#' This function identifies the point of maximum curvature (the "knee")
#' in a sorted signal vector, often used as a threshold to distinguish
#' high signal from background noise.
#'
#' @param signal A numeric vector of signal values.
#' @param plot Logical, if TRUE a diagnostic plot will be shown with the detected knee.
#'
#' @return A list with:
#'   \item{knee_index}{The index of the knee in the sorted signal vector}
#'   \item{threshold_value}{The signal value at the knee (used as threshold)}
#'
#' @examples
#' result <- find_signal_knee(runif(1000), plot = TRUE)
#' result$threshold_value
find_signal_knee <- function(signal, plot = FALSE) {
  # Sort unique signal values (for knee detection)
  signal <- sort(unique(signal))
  
  # Create x and y axes
  x <- seq_along(signal)
  y <- signal
  
  # Define a line from the first to last point
  point1 <- c(x[1], y[1])
  point2 <- c(x[length(x)], y[length(y)])
  
  # Normalize the vector representing the line
  line_vec <- point2 - point1
  line_vec_norm <- line_vec / sqrt(sum(line_vec^2))
  
  # Vectors from point1 to all other points
  vecs_to_line <- cbind(x - point1[1], y - point1[2])
  
  # Project each vector onto the line
  proj_lengths <- rowSums(vecs_to_line * matrix(rep(line_vec_norm, length(x)), ncol = 2, byrow = TRUE))
  proj_points <- matrix(rep(point1, length(x)), ncol = 2, byrow = TRUE) +
    proj_lengths %*% t(line_vec_norm)
  
  # Calculate perpendicular distances from each point to the line
  distances <- sqrt(rowSums((cbind(x, y) - proj_points)^2))
  
  # Find the index with maximum distance ? this is the "knee"
  knee_index <- which.max(distances)
  threshold_value <- signal[knee_index]
  
  # Optional plot to visualize the knee
  if (plot) {
    plot(x, y, type = "l", main = "Knee Detection", xlab = "Index", ylab = "Signal")
    points(x[knee_index], y[knee_index], col = "red", pch = 19)
  }
  
  # Return both the index and the signal threshold
  list(knee_index = knee_index, threshold_value = threshold_value)
}

