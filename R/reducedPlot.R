#' @title Plotting function for reduced dimensions
#' @description This function generates a plot for reduced dimensions of single-cell RNA-seq data
#' @param sce A SingleCellExperiment object containing the single-cell RNA-seq data.
#' @param var The variable to plot on the reduced dimensions.
#' @param axis_title_size Size of the axis titles (default is 13).
#' @param axis_text_size Size of the axis text (default is 11).
#' @param legend_text_size Size of the legend text (default is 11).
#' @param legend_title_size Size of the legend title (default is 13).
#' @param title_size Size of the plot title (default is 15).
#' @param ... Additional arguments passed to the dittoDimPlot function.
#' @return A ggplot object representing the reduced dimension plot.
#' @importFrom ggplot2 theme_classic theme element_text
#' @export
reducedPlot <- function(sce, var, axis_title_size = 15, axis_text_size = 14,
                        legend_text_size = 14, legend_title_size = 15,
                        title_size = 15, ...) {
  dittoSeq::dittoDimPlot(sce, var, ...) + theme_classic() +
    theme(
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size),
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_title_size),
      title = element_text(size = title_size)
    )
}



