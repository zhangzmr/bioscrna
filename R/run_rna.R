#' @title Process RNA data
#' @description This function processes RNA data in a SingleCellExperiment object.
#' It normalizes the data, performs dimension reduction, and merges batches using specified methods.
#' @param sce A SingleCellExperiment object containing the RNA data.
#' @param batch A character string specifying the batch column in the SingleCellExperiment object (default is "batch").
#' @param dim An integer specifying the number of dimensions to use for correction (default is 50).
#' @param n An integer specifying the number of high variable genes to select for PCA (default is 2000).
#' @param rank An integer specifying the number of principal components to compute (default is 50).
#' @param filter_genes A character vector of gene names to filter out from the analysis (default is c("IGLV", "IGHV", "IGKV", "IGLC")).
#' @param merge_method A character vector specifying the method to use for merging batches (default is c("FastMNN", "Harmony")).
#' @param cores An integer specifying the number of cores to use for parallel processing (default is 8).
#' @return A SingleCellExperiment object with processed RNA data, including normalized counts, dimension reduction, and batch correction.
#' @export


### Process RNA data

run_rna <- function(sce, batch = "batch", dim = 50, n = 2000, rank = 50,
                    filter_genes = c("IGLV", "IGHV", "IGKV", "IGLC"),
                    merge_method = c("FastMNN", "Harmony"), cores = 8) {
  sce <- scranNormalize(sce, graph.fun = "louvain")
  if (is.null(filter_genes)) {
    print("Use all gene for downstream analysis")
    temp <- sce
  } else {
    index <- c()
    for (i in filter_genes) {
      ind <- grep(i, rownames(sce))
      index <- c(index, ind)
      print(i)
    }
    genes <- rownames(sce)[index]
    print(paste0("Filtering out ", length(genes), " genes"))
    temp <- sce[!rownames(sce) %in% genes, ]
  }
  print(paste0("Running dimension reduction on ", ncol(temp), " cells and ", nrow(temp), " genes"))
  temp <- dimensionReduction(temp, n = n, rank = rank)
  temp <- merging(temp, batch = batch, dim = dim, method = merge_method, cores = cores)

  ## add dimension reduction data
  print("Adding dimension reduction data to the SCE object")
  # check if PCA and UMAP are present in temp
  if (!"PCA" %in% reducedDimNames(temp)) {
    stop("PCA is not present in the temp object")
  } else {
    reducedDim(sce, "PCA") <- reducedDim(temp, "PCA")
  }
  if (!"UMAP" %in% reducedDimNames(temp)) {
    stop("UMAP is not present in the temp object")
  } else {
    reducedDim(sce, "UMAP") <- reducedDim(temp, "UMAP")
  }

  if (!"batch_corrected" %in% reducedDimNames(temp)) {
    message("batch_corrected is not present in the temp object")
  } else {
    reducedDim(sce, "batch_corrected") <- reducedDim(temp, "batch_corrected")
  }
  if (!"UMAP_corrected" %in% reducedDimNames(temp)) {
    message("UMAP_corrected is not present in the temp object")
  } else {
    reducedDim(sce, "UMAP_corrected") <- reducedDim(temp, "UMAP_corrected")
  }

  if (!"HARMONY" %in% reducedDimNames(temp)) {
    message("HARMONY is not present in the temp object")
  } else {
    reducedDim(sce, "HARMONY") <- reducedDim(temp, "HARMONY")
  }
  if (!"UMAP_HARMONY" %in% reducedDimNames(temp)) {
    message("UMAP_HARMONY is not present in the temp object")
  } else {
    reducedDim(sce, "UMAP_HARMONY") <- reducedDim(temp, "UMAP_HARMONY")
  }
  return(sce)
}
