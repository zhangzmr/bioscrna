#' @title Clustering function for SingleCellExperiment objects
#' @param sce A SingleCellExperiment object.
#' @param k_range A vector of integers specifying the range of k values to test for clustering.
#' @param min_cl Minimum number of clusters to consider.
#' @param method Clustering method to use, e.g., "louvain", "kmeans", etc.
#' @param dimred Dimension reduction method to use, e.g., "PCA", "UMAP".
#' @param cores Number of cores to use for parallel processing.
#' @param top Number of top clusters to return based on silhouette score.
#' @return A list with \code{sce} (SingleCellExperiment with clustering results in colData)
#'   and \code{metrics} (data.frame with num.clusters, silhouette, and wcss per parameter set).
#' @importFrom BiocParallel MulticoreParam
#' @importFrom bluster NNGraphParam
#' @export
clustering <- function(sce, k_range, min_cl = 10, method = "louvain", dimred = "PCA", cores = 8, top = 5) {
  message(paste0("Clustering method using: ", method))
  message(paste0("Reduced Dimension use: ", dimred))

  print("Clustering with different K values...")
  out <- bluster::clusterSweep(reducedDim(sce, dimred), NNGraphParam(),
                      k=as.integer(k_range),
                      cluster.fun=method,
                      BPPARAM=BiocParallel::MulticoreParam(cores))

  message(paste0("Finding the best top ", top, " clusters..."))
  df <- as.data.frame(out$parameters)
  df$num.clusters <- vapply(as.list(out$clusters), function(cluster) {
    length(unique(cluster))
  }, 0L)

  all.sil <- lapply(as.list(out$clusters), function(cluster) {
    sil <- bluster::approxSilhouette(reducedDim(sce, dimred), cluster)
    mean(sil$width)
  })
  df$silhouette <- unlist(all.sil)

  all.wcss <- lapply(as.list(out$clusters), function(cluster) {
    sum(bluster::clusterRMSD(reducedDim(sce, dimred), cluster, sum = TRUE), na.rm = TRUE)
  })
  df$wcss <- unlist(all.wcss)

  df <- df[df$num.clusters >= min_cl, ]
  df <- df[order(-df$silhouette), ]

  ## get clusters
  clname <- rownames(df[1:top, ])
  clustering <- data.frame(out$clusters)
  clustering <- clustering[ , clname]

  name <- grep("clust", colnames(colData(sce)))
  if (length(name) > 0) {
    message("Replacing old clustering results...")
    colData(sce) <- colData(sce)[ , -name]
  }
  colData(sce) <- cbind(colData(sce), clustering)
  return(list(sce = sce, metrics = df))
}


