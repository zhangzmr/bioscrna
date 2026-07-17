#' @title Merging Function
#' @description This function merges data from a SingleCellExperiment object using specified methods.
#' It supports FastMNN and Harmony methods for batch correction and dimensionality reduction.
#' It requires the `batchelor`, `harmony`, and `scater` packages for the analysis.
#' @param sce A SingleCellExperiment object containing the data.
#' @param batch A character string specifying the batch column in the SingleCellExperiment object.
#' @param dim An integer specifying the number of dimensions to use for correction.
#' @param method A character vector specifying the method to use for merging. Options are "FastMNN" and "Harmony".
#' @param cores An integer specifying the number of cores to use for parallel processing.
#' @return A SingleCellExperiment object after merging the data.
#' @importFrom BiocParallel MulticoreParam
#' @export
merging <- function(sce, batch = "batch", dim = 50,
                    method = c("FastMNN", "Harmony"), cores = 9) {
  if ("FastMNN" %in% method) {
    if (!requireNamespace("batchelor", quietly = TRUE)) {
      stop("Please install batchelor package.")
    }
    message("Running FastMNN...")
    sce_corrected <- batchelor::fastMNN(
      sce, batch = sce[[batch]], d = dim,
      BPPARAM = MulticoreParam(cores)
    )
    reducedDim(sce, "batch_corrected") <- reducedDim(sce_corrected, "corrected")
    sce <- scater::runUMAP(
      sce, dimred = "batch_corrected", name = "UMAP_corrected",
      BPPARAM = MulticoreParam(cores)
    )
  }

  if ("Harmony" %in% method) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("Please install the harmony package to use method = \"Harmony\".")
    }

    message("Running Harmony...")
    sce <- harmony::RunHarmony(sce, batch, ncores = cores, verbose = TRUE)
    sce <- scater::runUMAP(
      sce, dimred = "HARMONY", name = "UMAP_HARMONY",
      BPPARAM = MulticoreParam(cores)
    )
  }
  return(sce)
}


