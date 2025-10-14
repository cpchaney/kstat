library(Matrix)
library(imager)
library(MatrixGenerics)

processLandmarks <- function(readsFile, width, height, flipped.horizontal = FALSE, scale.factor = 1) {
  # Read the reads file
  reads <- read.csv(file = readsFile)
  landmarks <- unique(reads$gene)

  # Scale and flip coordinates as needed
  reads$X <- floor(reads$X / scale.factor)
  if (flipped.horizontal) {
    reads$X <- width - reads$X + 1
  }
  reads$Y <- height - floor(reads$Y / scale.factor)

  # Create sparse matrices for each landmark
  landmarks_matrix <- sapply(landmarks, function(landmark) {
    landmark_reads <- subset(reads, gene == landmark)
    landmark_matrix <- sparseMatrix(
      i = landmark_reads$Y,
      j = landmark_reads$X,
      dims = c(height, width)
    )
    as.vector(t(landmark_matrix), mode = "integer")
  })

  # Assign row names and return as `dgCMatrix`
  row.names(landmarks_matrix) <- landmarks
  return(as(landmarks_matrix, "dgCMatrix"))
}

getMask <- function(mask_file, height, width, invert = FALSE) {
  # Load the mask image
  mask_image <- load.image(mask_file)
  if (invert) {
    mask_image <- 1 - mask_image
  }

  # Convert to sparse matrix
  mask <- as(t(as.matrix(mask_image)), "dgTMatrix")
  mask@x <- pmax(0, pmin(1, mask@x)) # Binarize mask values

  # Transform coordinates
  return(sparseMatrix(i = height - mask@i, j = mask@j + 1, x = mask@x, dims = c(height, width)))
}

getExpressionData <- function(sceFile, assay_name = "imputed") {
  # Load single-cell experiment data
  return(assay(readRDS(file = sceFile), assay_name))
}

encodeBin <- function(bin) {
  # Encode binary vector as a number
  return(sum(bin * 2^(seq_along(bin) - 1)))
}

getBinCoordinates <- function(bin, height, width) {
  # Calculate x-y coordinates from bin
  j <- bin %% width
  i <- bin %/% width + 1
  j[j == 0] <- width
  return(cbind(i, j))
}

getSignalCoordinates <- function(channel_name, v, height, width, mask = NULL) {
  # Apply mask if provided
  if (!is.null(mask)) {
    v <- v * mask
  }

  # Create sparse matrix from vector
  M <- as(matrix(v, nrow = height, ncol = width, byrow = TRUE), "TsparseMatrix")

  # Create data frame with coordinates and signal magnitudes
  return(data.frame(
    channel = factor(rep(channel_name, length(M@x))),
    X = M@j + 1,
    Y = height - M@i + 1,
    magnitude = M@x
  ))
}

getLandmarkSignal <- function(landmark, readsFile, height) {
  # Filter reads for the specified landmark
  reads <- read.csv(file = readsFile)
  reads$X <- floor(reads$X)
  reads$Y <- height - floor(reads$Y)
  return(subset(reads, gene == landmark))
}

getSpatialSignal <- function(signal, cell_probabilities, mapped.bins = NULL, log = TRUE) {
  # Calculate weighted signal
  if (log) {
    weighted_signal <- log(((exp(signal) - 1) %*% cell_probabilities) + 1)
  } else {
    weighted_signal <- signal %*% cell_probabilities
  }

  # Adjust for mapped bins if provided
  if (!is.null(mapped.bins)) {
    weighted_signal <- weighted_signal[, mapped.bins, drop = FALSE]
  }

  # Normalize signal
  if (is.null(dim(signal))) {
    weighted_signal <- weighted_signal * max(signal) / max(weighted_signal)
  } else {
    weighted_signal <- weighted_signal * (rowMaxs(as.matrix(signal)) / rowMaxs(as.matrix(weighted_signal)))
  }

  return(weighted_signal)
}

gaussianBlur <- function(spatial_signal, height, width, sigma = 1) {
  ski <- import("skimage")
  # Convert signal into a matrix
  spatial_matrix <- matrix(spatial_signal, nrow = height, ncol = width, byrow = TRUE)

  # Apply Gaussian blur
  blurred <- ski$filters$gaussian(spatial_matrix, sigma = sigma)

  # Normalize and return
  ## row_maximums <- rowMaxs(blurred)
  ## row_maximums[row_maximums == 0] <- 1
  ## blurred <- blurred * (rowMaxs(as.matrix(spatial_signal)) / row_maximums)
  return(as.numeric(t(blurred)))
}

splitMatrixIntoTiles <- function(matrix, tile_size) {
  # Split a matrix into tiles of specified size
  M <- nrow(matrix)
  N <- ncol(matrix)
  num_tiles_row <- ceiling(M / tile_size)
  num_tiles_col <- ceiling(N / tile_size)

  tiles <- list()
  for (i in seq_len(num_tiles_row)) {
    for (j in seq_len(num_tiles_col)) {
      row_start <- (i - 1) * tile_size + 1
      row_end <- min(i * tile_size, M)
      col_start <- (j - 1) * tile_size + 1
      col_end <- min(j * tile_size, N)
      tiles[[paste(i, j, sep = "_")]] <- matrix[row_start:row_end, col_start:col_end]
    }
  }
  return(tiles)
}

collapseBins <- function(signal, height, width, tile_size = 4) {
  # Reduce resolution by collapsing bins
  signal_matrix <- matrix(signal, nrow = height, ncol = width, byrow = TRUE)
  horizontal_padding <- tile_size - (width %% tile_size)
  signal_matrix <- cbind(signal_matrix, matrix(0, nrow = height, ncol = horizontal_padding))
  vertical_padding <- tile_size - (height %% tile_size)
  signal_matrix <- rbind(signal_matrix, matrix(0, nrow = vertical_padding, ncol = ncol(signal_matrix)))

  tiles <- splitMatrixIntoTiles(signal_matrix, tile_size)
  return(list(
    signal = sapply(tiles, sum, simplify = TRUE),
    height = nrow(signal_matrix) / tile_size,
    width = ncol(signal_matrix) / tile_size
  ))
}
