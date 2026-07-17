#' @title scranNormalize
#' @description Normalize single-cell RNA-seq data using scran package.
#' @param sce A SingleCellExperiment object containing the data.
#' @param ... Additional arguments passed to the `scran::quickCluster` function.
#' @return A SingleCellExperiment object after normalization.
#' @importFrom BiocParallel MulticoreParam
#' @export
scranNormalize <- function(sce, ...) {
  message("Run quick cluster...")
  clust <- scran::quickCluster(sce, ...)
  message("Run computer sum factors...")
  sce <- scran::computeSumFactors(sce, clusters = clust, BPPARAM = MulticoreParam())
  message("Run log norm...")
  sce <- scuttle::logNormCounts(sce)
  return(sce)
}

