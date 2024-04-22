#' Perform Principal Component Analysis (PCA)
#'
#' This function performs Principal Component Analysis (PCA) on the input data matrix.
#'
#' @importFrom stats prcomp
#' @param data Input data matrix.
#'
#' @return PCA result object.
#'
#' @export
perform_pca <- function(data) {
  # Extract expression data
  expression_data <- data[,1:ncol(data)]

  # Perform PCA
  pca_result <- prcomp(t(expression_data))
  return(pca_result)
}
