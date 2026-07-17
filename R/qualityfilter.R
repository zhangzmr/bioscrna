#' @title Quality Control for Single-Cell RNA-seq Data
#' @description This function performs quality control on a SingleCellExperiment object.
#' It filters cells based on mitochondrial gene percentage, total counts, number of detected genes,
#' and optionally filters out lowly expressed genes.
#' It also allows for the removal of specific quality control genes such as mitochondrial and ribosomal genes.
#' @param sce A SingleCellExperiment object containing the data.
#' @param mito_percent Maximum percentage of mitochondrial genes allowed in a cell.
#' @param filter_outlier Whether to filter out outlier cells based on total counts.
#' @param sum_max Maximum total counts allowed for a cell.
#' @param sum_min Minimum total counts required for a cell.
#' @param gene_max Maximum number of detected genes allowed for a cell.
#' @param gene_min Minimum number of detected genes required for a cell.
#' @param filter_gene Whether to filter out lowly expressed genes.
#' @param min_counts Minimum counts required for a gene to be retained.
#' @param ... Additional arguments passed to the `scuttle::isOutlier` function
#' @return A SingleCellExperiment object after quality control.
#' @importFrom MatrixGenerics rowSums2
#' @export
qualityfilter <- function(sce, mito_percent = 10, filter_outlier = TRUE,
                          sum_max = 16000, sum_min = 1000, gene_max = 4000,
                          gene_min = 700, filter_gene = TRUE, min_counts = 5, ...) {
  if (isTRUE(filter_outlier)) {
    message("Find outliers...")
    sum.outlier <- scuttle::isOutlier(sce$sum, type = "both", log = TRUE, ...)
    print(table(sum.outlier))
    sce <- sce[, !sum.outlier]
  }

  message("Filter with fix threshold...")
  sum.max <- sce$sum < sum_max
  sum.min <- sce$sum > sum_min
  detected.max <- sce$detected < gene_max
  detected.min <- sce$detected > gene_min
  mito_cell <- sce$subsets_Mito_percent < mito_percent
  keep <- sum.max & sum.min & detected.max & detected.min & mito_cell
  print(table(keep))
  sce <- sce[, keep]

  print("Filtering mitochondrial genes, ribosomal genes...")
  sce <- removeQCGenes(sce)

  if (isTRUE(filter_gene)) {
    message("Filter lowly expressed genes...")
    rowcounts <- rowSums2(assay(sce, "counts"))
    keep_genes <- names(rowcounts)[rowcounts >= min_counts]
    sce <- sce[rownames(sce) %in% keep_genes, ]
  } else {
    message("No lowly expressed gene filtering applied.")
  }
  return(sce)
}

#' @title Remove common QC genes
#' @description Remove mitochondrial, ribosomal, and MALAT1/Malat1 genes from an SCE.
#' @param sce A SingleCellExperiment object.
#' @param type Organism type, either `"Human"` or `"Mouse"`.
#' @return A SingleCellExperiment object with QC genes removed.
#' @export
removeQCGenes <- function(sce, type = "Human") {
  if (type == "Human") {
    mito <- rownames(sce)[grep("MT-", rownames(sce))]
    ribo <- rownames(sce)[grep("^RPL|^RPS", rownames(sce))]
    remove_genes <- c("MALAT1", mito, ribo)
  } else if (type == "Mouse") {
    mito <- rownames(sce)[grep("mt-", rownames(sce))]
    ribo <- rownames(sce)[grep("^Rpl|^Rps", rownames(sce))]
    remove_genes <- c("Malat1", mito, ribo)
  } else {
    stop("type must be \"Human\" or \"Mouse\".")
  }

  sce <- sce[!rownames(sce) %in% remove_genes, ]
  return(sce)
}
