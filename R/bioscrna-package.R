#' @keywords internal
"_PACKAGE"

#' @importFrom utils read.delim
#' @importFrom stats dist hclust
#' @importFrom methods as
#' @importFrom grDevices pdf dev.off
#' @importFrom SummarizedExperiment assay assayNames colData rowData
#' @importFrom SummarizedExperiment assay<- colData<- rowData<-
#' @importFrom SingleCellExperiment SingleCellExperiment counts reducedDim reducedDimNames
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom MatrixGenerics rowSums2
#' @importFrom BiocParallel MulticoreParam
#' @importFrom bluster NNGraphParam
NULL

utils::globalVariables(c("Group", "Feature", "PctExpr", "AvgExpr_capped"))
