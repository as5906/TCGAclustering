#' Visualize PCA results
#'
#' This function visualizes the results of Principal Component Analysis (PCA) using a scatter plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text labs ggtitle
#' @param pca_result PCA result object.
#' @param color_data Optional vector specifying color for data points.
#' @param labels Optional vector of labels for data points.
#'
#' @return A ggplot object representing the PCA scatter plot.
#'
#' @export
graph_pca <- function(pca_result, color_data=NULL, labels=NULL) {
  # Extract PC1 and PC2 scores
  pc1 <- pca_result$x[,1]
  pc2 <- pca_result$x[,2]

  # Create a data frame with PC1 and PC2 scores
  pca_df <- data.frame(PC1 = pc1, PC2 = pc2)

  # Assign color if provided
  if (!is.null(color_data)) {
    pca_df$Color <- color_data
  }

  # If color_data is not provided, assign a default color to pca_df$Color
  else {
    pca_df$Color <- "black"
  }

  # Plot PC1 vs PC2
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Color))
  p <- p + geom_point()
  p <- p + labs(x = "PC1", y = "PC2")

  # Add labels if provided
  if (!is.null(labels)) {
    p <- p + geom_text(label = labels, vjust = -0.5)
  }

  return(p)
}
