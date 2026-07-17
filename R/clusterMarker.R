#' @title Cluster Marker Detection
#' @description This function identifies top markers for each cluster in a SingleCellExperiment object.
#' It uses the `scran` package to find markers based on the specified clustering variable.
#' The function returns a data frame containing the top markers for each cluster, filtered by FDR and sorted by log fold change.
#' @param sce A SingleCellExperiment object.
#' @param top Number of top markers to return for each cluster.
#' @param cluster The name of the clustering variable in colData(sce).
#' @return A data frame containing the top markers for each cluster.
#' @importFrom BiocParallel MulticoreParam
#' @export
clusterMarker <- function(sce, top = 5, cluster) {
  marker <- scran::findMarkers(sce, groups = cluster, pval.type = "some",
                               direction = "up", BPPARAM = MulticoreParam())
  markerTop <- data.frame()
  for (l in names(marker@listData)) {
    temp <- as.data.frame(marker@listData[[l]])
    temp <- temp[temp$FDR < 0.05,]
    temp <- temp[order(-temp$summary.logFC), ]
    temp <- temp[1:top, 1:3]
    temp$cluster <- l
    temp$genes <- rownames(temp)
    rownames(temp) <- NULL
    markerTop <- rbind(markerTop, temp)
  }
  return(markerTop)
}