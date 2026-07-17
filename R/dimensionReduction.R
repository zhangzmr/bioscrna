#' @title Dimension Reduction
#' @description This function performs dimension reduction on a SingleCellExperiment object using PCA and UMAP.
#' It selects high variable genes, computes PCA, and then runs UMAP for visualization.
#' It requires the `scran` and `scater` packages for the analysis.
#' @param sce A SingleCellExperiment object containing the data.
#' @param n Number of high variable genes to select for PCA.
#' @param rank Number of principal components to compute.
#' @param dimred Type of dimension reduction to perform (default is "PCA").
#' @param ... Additional arguments passed to the `scater::runUMAP` function
#' @return A SingleCellExperiment object with reduced dimensions.
#' @importFrom BiocParallel MulticoreParam
#' @export
dimensionReduction <- function(sce, n = 1000, rank = 50, dimred = "PCA", ...) {
  message("Choosing High Variable Genes...")
  dec <- scran::modelGeneVar(sce)
  hvg <- scran::getTopHVGs(dec, n = n)
  message("Running PCA...")
  sce <- scran::fixedPCA(sce, rank = rank, subset.row = hvg, assay.type = "logcounts")
  message("Running UMAP...")
  sce <- scater::runUMAP(sce, dimred = dimred, BPPARAM = MulticoreParam(), ...)
  return(sce)
}

