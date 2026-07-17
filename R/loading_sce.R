#' @title Combine Multiple SingleCellExperiment Objects
#' @description Combines a list of SingleCellExperiment objects by taking the intersection of genes
#' and merging count matrices, column data, and row data.
#' @param sce_list A named list of SingleCellExperiment objects to combine.
#' @param sample A character vector of sample names corresponding to the SCE objects.
#' @return A combined SingleCellExperiment object with all samples.
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @keywords internal
combine_sce_list <- function(sce_list, sample) {
  message("Combining gene expressions...")
  
  # Get common genes across all samples (intersection)
  all_genes <- lapply(sce_list, rownames)
  common_genes <- Reduce(intersect, all_genes)
  message(paste0("Found ", length(common_genes), " common genes across all samples"))
  
  # Subset each SCE to common genes and combine count matrices
  combined_counts_list <- list()
  combined_coldata_list <- list()
  combined_rowdata_list <- list()
  
  for (i in 1:length(sce_list)) {
    sce_temp <- sce_list[[i]]
    # Subset to common genes
    sce_temp <- sce_temp[common_genes, ]
    # Extract counts matrix
    combined_counts_list[[i]] <- counts(sce_temp)
    # Extract and modify colData
    temp_coldata <- colData(sce_temp)
    temp_coldata$batch <- sample[i]
    combined_coldata_list[[i]] <- temp_coldata
    # Store rowData (should be the same for common genes)
    if (i == 1) {
      combined_rowdata_list <- rowData(sce_temp)
    }
  }
  
  # Combine count matrices by columns (cbind)
  combined_counts <- do.call(cbind, combined_counts_list)
  combined_coldata <- do.call(rbind, combined_coldata_list)
  
  # Create combined SCE object
  sce <- SingleCellExperiment(assays = list(counts = combined_counts),
                             colData = combined_coldata,
                             rowData = combined_rowdata_list)
  
  message("Done combining...")
  return(sce)
}

#' @title Loading SingleCellExperiment Object
#' @description This function loads 10X Genomics data into a SingleCellExperiment object.
#' It reads the data from specified directories and combines them if multiple directories are provided.
#' It also adds quality control metrics for mitochondrial and ribosomal genes.
#' @details The function reads 10X data manually and combines SingleCellExperiment objects by taking the intersection of genes across samples.
#' It checks if the length of the `directory` and `sample` arguments are the same, stopping if they are not.
#' @param directory A character vector of directories containing 10X Genomics data.
#' @param sample A character vector of sample names corresponding to the directories.
#' @param add_qc Logical, whether to add quality control metrics for mitochondrial and ribosomal genes. Default is TRUE.
#' @return A SingleCellExperiment object containing the loaded data.
#' @export

loading_sce <- function(directory, sample, add_qc = TRUE) {
  if (length(directory) != length(sample)) {
    stop("The length of directory and sample names should be the same.")
  }

  if (length(directory) == 1) {
    message("Number of samples is 1, reading data...")
    sce <- loading_oneSCE(directory, sample)
    sce$batch <- sce$Sample
  }

  if (length(directory) > 1) {
    message(paste0("Number of samples is ", length(directory), ", reading data..."))
    sce_list <- list()
    for (i in 1:length(directory)) {
      dir <- directory[i]
      name <- sample[i]
      temp <- loading_oneSCE(dir, name)
      sce_list[[name]] <- temp
      message(paste0("Finish reading data: ", name))
    }
    
    # Combine all SCE objects
    sce <- combine_sce_list(sce_list, sample)
  }

  # add qc stats if requested
  if (add_qc) {
    message("Adding quality control metrics...")
    is.mito <- grep("MT-", rownames(sce))
    is.ribo <- grep("^RPL|^RPS", rownames(sce))
    sce <- scuttle::addPerCellQCMetrics(sce, subsets=list(Mito=is.mito, Ribo=is.ribo))
  }
  return(sce)
}



