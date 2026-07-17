#' Clustered Dot Plot for SingleCellExperiment
#' This function creates a clustered dot plot for visualizing gene expression
#' across different groups in a SingleCellExperiment object.
#' @param sce A SingleCellExperiment object containing the data.
#' @param features A character vector of gene names to include in the plot.
#' @param group_by A string specifying the column in colData(sce) to group cells by.
#' @param assay_name The name of the assay to use for expression data (default is "logcounts").
#' @param cluster_features Logical, whether to cluster features (genes) (default is TRUE).
#' @param cluster_groups Logical, whether to cluster groups (default is TRUE).
#' @param scale_expression Logical, whether to scale expression values (default is TRUE).
#' @param dot_min Minimum percentage for dot size scaling (default is 0).
#' @param dot_max Maximum percentage for dot size scaling (default is 100).
#' @param col_min Minimum expression value for color scaling (default is -2.5).
#' @param col_max Maximum expression value for color scaling (default is 2.5).
#' @param colors A vector of three colors for low, mid, and high expression (default is c("blue", "white", "red")).
#' @param font_size Font size for axis text and titles (default is 10).
#' @param angle Angle for x-axis text (default is 45).
#' @param title Optional title for the plot.
#' @return A ggplot2 object representing the clustered dot plot.
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous scale_color_gradient2 theme_classic theme element_text element_blank guides guide_legend guide_colorbar ggtitle
#' @export
clustered_dotplot <- function(sce, features, group_by,
                              assay_name = "logcounts",
                              cluster_features = TRUE,
                              cluster_groups = TRUE,
                              scale_expression = TRUE,
                              dot_min = 0, dot_max = 100,
                              col_min = -2.5, col_max = 2.5,
                              colors = c("blue", "white", "red"),
                              font_size = 10,
                              angle = 45,
                              title = NULL) {
  
  # Validate inputs
  if (!group_by %in% colnames(colData(sce))) {
    stop(paste("Group variable", group_by, "not found in colData"))
  }
  
  # Filter features that exist in the data
  features_available <- intersect(features, rownames(sce))
  if (length(features_available) == 0) {
    stop("None of the specified features found in the data")
  }
  
  if (length(features_available) < length(features)) {
    missing_features <- setdiff(features, features_available)
    warning(paste("Features not found:", paste(missing_features, collapse = ", ")))
  }
  
  # Get expression data
  if (assay_name %in% assayNames(sce)) {
    expr_data <- assay(sce, assay_name)
  } else {
    # If logcounts doesn't exist, use counts and log-transform
    warning("logcounts not found, using log-transformed counts")
    expr_data <- log1p(assay(sce, "counts"))
  }
  
  # Calculate average expression and percentage expressed for each group
  groups <- factor(sce[[group_by]])
  group_levels <- levels(groups)
  
  # Initialize matrices
  avg_expr <- matrix(0, nrow = length(features_available), ncol = length(group_levels))
  pct_expr <- matrix(0, nrow = length(features_available), ncol = length(group_levels))
  
  rownames(avg_expr) <- features_available
  colnames(avg_expr) <- group_levels
  rownames(pct_expr) <- features_available
  colnames(pct_expr) <- group_levels
  
  # Calculate statistics for each group
  for (i in seq_along(group_levels)) {
    group_cells <- which(groups == group_levels[i])
    if (length(group_cells) > 0) {
      # Extract expression data for this group - ALWAYS use drop = FALSE
      group_expr <- expr_data[features_available, group_cells, drop = FALSE]
      
      # Ensure group_expr is always a matrix
      if (!is.matrix(group_expr)) {
        group_expr <- as.matrix(group_expr)
      }
      
      # Handle single gene case specially
      if (length(features_available) == 1) {
        # Single gene case
        if (ncol(group_expr) == 1) {
          # Single gene, single cell
          avg_expr[, i] <- group_expr[1, 1]
          pct_expr[, i] <- (group_expr[1, 1] > 0) * 100
        } else {
          # Single gene, multiple cells
          avg_expr[, i] <- mean(group_expr[1, ])
          pct_expr[, i] <- sum(group_expr[1, ] > 0) / ncol(group_expr) * 100
        }
      } else {
        # Multiple genes case
        if (ncol(group_expr) == 1) {
          # Multiple genes, single cell
          avg_expr[, i] <- group_expr[, 1]
          pct_expr[, i] <- (group_expr[, 1] > 0) * 100
        } else {
          # Multiple genes, multiple cells
          avg_expr[, i] <- rowMeans(group_expr)
          pct_expr[, i] <- rowSums(group_expr > 0) / ncol(group_expr) * 100
        }
      }
    }
  }
  
  # Scale expression if requested
  if (scale_expression) {
    avg_expr <- t(scale(t(avg_expr)))
    avg_expr[is.na(avg_expr)] <- 0
  }
  
  # Cluster features (genes) if requested
  if (cluster_features && length(features_available) > 1) {
    feature_dist <- dist(avg_expr)
    feature_clust <- hclust(feature_dist, method = "ward.D2")
    features_available <- features_available[feature_clust$order]
  }
  
  # Cluster groups if requested
  if (cluster_groups && length(group_levels) > 1) {
    group_dist <- dist(t(avg_expr))
    group_clust <- hclust(group_dist, method = "ward.D2")
    group_levels <- group_levels[group_clust$order]
  }
  
  # Prepare data for plotting
  plot_data <- expand.grid(Feature = features_available, Group = group_levels, stringsAsFactors = FALSE)
  
  # Extract values properly by matching feature and group names
  plot_data$AvgExpr <- numeric(nrow(plot_data))
  plot_data$PctExpr <- numeric(nrow(plot_data))
  
  for (i in 1:nrow(plot_data)) {
    feat <- as.character(plot_data$Feature[i])
    grp <- as.character(plot_data$Group[i])
    plot_data$AvgExpr[i] <- avg_expr[feat, grp]
    plot_data$PctExpr[i] <- pct_expr[feat, grp]
  }
  
  # Cap percentage values
  plot_data$PctExpr <- pmax(pmin(plot_data$PctExpr, dot_max), dot_min)
  
  # Cap expression values for color scale
  plot_data$AvgExpr_capped <- pmax(pmin(plot_data$AvgExpr, col_max), col_min)
  
  # Set factor levels to maintain order
  plot_data$Feature <- factor(plot_data$Feature, levels = features_available)
  plot_data$Group <- factor(plot_data$Group, levels = group_levels)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Group, y = Feature)) +
    geom_point(aes(size = PctExpr, color = AvgExpr_capped)) +
    scale_size_continuous(name = "Percent\nExpressed",
                          range = c(0, 6),
                          limits = c(dot_min, dot_max),
                          breaks = c(0, 25, 50, 75, 100)) +
    scale_color_gradient2(name = if(scale_expression) "Avg Expr\n(Scaled)" else "Avg Expr",
                          low = colors[1], mid = colors[2], high = colors[3],
                          midpoint = 0,
                          limits = c(col_min, col_max)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1, size = font_size),
      axis.text.y = element_text(size = font_size),
      axis.title = element_blank(),
      legend.title = element_text(size = font_size),
      legend.text = element_text(size = font_size - 1),
      plot.title = element_text(size = font_size + 2, hjust = 0.5)
    ) +
    guides(
      size = guide_legend(override.aes = list(color = "black")),
      color = guide_colorbar(order = 2)
    )
  
  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}


