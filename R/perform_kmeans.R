#' Perform K-means clustering with silhouette analysis
#'
#' This function performs K-means clustering with silhouette analysis to determine the optimal number of clusters.
#'
#' @importFrom cluster silhouette
#' @importFrom graphics abline
#' @importFrom stats complete.cases dist kmeans prcomp
#' @param pca_result PCA result object.
#' @param max_clusters Maximum number of clusters to consider.
#'
#' @return K-means clustering result object.
#'
#' @export
perform_kmeans_with_silhouette <- function(pca_result, max_clusters) {
  # Extract principal components
  pcs <- pca_result$x

  # Initialize vector to store silhouette widths
  silhouette_widths <- vector(mode = "numeric", length = max_clusters - 1)

  # Calculate silhouette widths for different numbers of clusters
  for (k in 2:max_clusters) {
    kmeans_result <- kmeans(pcs, centers = k)
    cluster_assignment <- kmeans_result$cluster
    sil_obj <- silhouette(cluster_assignment, dist(pcs))
    sil_width <- sil_obj[, "sil_width"]
    silhouette_widths[k - 1] <- mean(sil_width)
  }

  # Find the optimal K based on silhouette analysis
  optimal_k_silhouette <- which.max(silhouette_widths) + 1
  cat("Optimal number of clusters (K) suggested by silhouette analysis:", optimal_k_silhouette, "\n")

  # Plot silhouette widths
  plot(2:max_clusters, silhouette_widths, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of Clusters", ylab = "Average Silhouette Width",
       main = "Silhouette Analysis for Optimal K")
  abline(v = optimal_k_silhouette, col = "red", lty = 2)  # Add a vertical line for optimal K

  # Perform k-means clustering with the optimal K
  kmeans_result <- kmeans(pcs, centers = optimal_k_silhouette)

  return(kmeans_result)
}
