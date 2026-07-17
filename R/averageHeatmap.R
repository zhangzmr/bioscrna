#' @title Average Heatmap
#' @description This function generates an average heatmap for specified genes across groups in a SingleCell Experiment (SCE) object.
#' @param sce A SingleCellExperiment object containing the single-cell RNA-seq data.
#' @param genes A character vector of gene names to include in the heatmap.
#' @param groups A factor or character vector indicating the grouping of cells.
#' @param ... Additional arguments passed to the dittoHeatmap function.
#' @return A heatmap object showing the average expression of specified genes across groups.
#' @export

averageHeatmap <- function(sce, genes, groups, ...) {
  intergene <- intersect(genes, rownames(sce))
  average <- scuttle::sumCountsAcrossCells(sce, groups, exprs_values = "logcounts", average=TRUE)
  dittoSeq::dittoHeatmap(average, intergene, ...)
}



