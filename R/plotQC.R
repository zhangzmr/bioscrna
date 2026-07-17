#' @title Plot Quality Control
#' @description This function generates quality control plots for single-cell RNA-seq data.
#' @param sce A SingleCellExperiment object containing the single-cell RNA-seq data.
#' @param file The output file name for the PDF containing the plots.
#' @param width The width of the output PDF file.
#' @param height The height of the output PDF file.
#' @return The modified SingleCellExperiment object with quality control plots.
#' @importFrom ggplot2 theme_classic
#' @importFrom grDevices pdf dev.off
#' @export
plotQC <- function(sce, file = "plot/Quality_Control_before.pdf",
                   width = 7, height = 5) {
  message("Ploting started...")
  pdf(file, width = width, height = height)
  print(dittoSeq::dittoScatterPlot(
    object = sce, x.var = "sum", y.var = "detected",
    color.var = "batch", do.raster = TRUE
  ) + theme_classic())
  print(dittoSeq::dittoScatterPlot(
    object = sce, x.var = "sum", y.var = "subsets_Mito_percent",
    color.var = "batch", do.raster = TRUE
  ) + theme_classic())
  print(dittoSeq::dittoPlot(
    sce, "sum", group.by = "batch",
    plots = c("jitter", "boxplot"), jitter.size = 0.5, do.raster = TRUE
  ))
  print(dittoSeq::dittoPlot(
    sce, "detected", group.by = "batch",
    plots = c("jitter", "boxplot"), jitter.size = 0.5, do.raster = TRUE
  ))
  dev.off()
  message("Plot finished...")
  return(sce)
}

