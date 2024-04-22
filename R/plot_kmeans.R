#' Visualize K-means clustering results
#'
#' This function visualizes the results of K-means clustering using a scatter plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text labs ggtitle
#' @param pca_result PCA result object.
#' @param kmeans_result K-means clustering result object.
#'
#' @return A ggplot object representing the K-means clustering results.
#'
#' @export
plot_kmeans_results <- function(pca_result, kmeans_result) {
  # Extract principal components
  pcs <- pca_result$x

  # Create a dataframe with PC1 and PC2 scores and cluster labels
  df <- data.frame(PC1 = pcs[,1], PC2 = pcs[,2], Cluster = as.factor(kmeans_result$cluster))

  # Plot the clusters
  p <- ggplot(df, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point() +
    labs(x = "PC1", y = "PC2", color = "Cluster") +
    ggtitle("K-means Clustering Results")

  return(p)
}
