#' @title Plot Clusters
#' @description This function generates cluster plots for single-cell RNA-seq data.
#' @param sce A SingleCellExperiment object containing the single-cell RNA-seq data.
#' @param dimred The dimension reduction method to use for plotting (default is "UMAP").
#' @param file The output file name for the PDF containing the plots.
#' @param width The width of the output PDF file.
#' @param height The height of the output PDF file.
#' @param label Whether to label the clusters in the plot (default is TRUE).
#' @param size The size of the points in the plot (default is 0.8).
#' @param legend.size The size of the legend icons (default is 10).
#' @param labels.size The size of the cluster labels (default is 5).
#' @param do.raster Whether to rasterize the plot (default is TRUE).
#' @param ... Additional arguments passed to \code{\link[dittoSeq]{dittoDimPlot}}.
#' @return The modified SingleCellExperiment object with cluster plots.
#' @importFrom ggplot2 theme_classic theme element_text
#' @importFrom grDevices pdf dev.off
#' @export
plotCluster <- function(sce, dimred = "UMAP", file = "plot.pdf", width = 9, height = 8,
                        label = TRUE, size = 0.8, legend.size = 10, labels.size = 5,
                        do.raster = TRUE, ...) {
  plot_name <- colnames(colData(sce))[grep("clust", colnames(colData(sce)))]
  pdf(file, width = width, height = height)
  for (name in plot_name) {
    gra <- dittoSeq::dittoDimPlot(
      sce, name, reduction.use = dimred, size = size,
      do.label = label, legend.size = legend.size, labels.size = labels.size,
      do.raster = do.raster, ...
    ) + theme_classic() +
      theme(legend.text = element_text(size = 15),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15))
    print(gra)
  }
  dev.off()
  return(sce)
}

