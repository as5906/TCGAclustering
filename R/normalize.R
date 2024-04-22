#' Normalize and compute Z-scores
#'
#' This function performs log2 normalization, quantile normalization, and computes Z-scores for the input data.
#' @importFrom preprocessCore normalize.quantiles
#' @param data Input data matrix.
#'
#' @return A matrix of Z-scores.
#'
#' @export
normalize_and_zscore <- function(data) {
  # Log2 normalize each column
  log2_normalized_data <- log2(data + 1) # Adding 1 to avoid log(0)
  log2_normalized_matrix <- as.matrix(log2_normalized_data)

  # Quantile normalize each column
  quantile_normalized_data <- normalize.quantiles(log2_normalized_matrix, keep.names=TRUE)

  # Convert to z-scores for each column
  zscore_data <- scale(quantile_normalized_data)

  return(zscore_data)
}
