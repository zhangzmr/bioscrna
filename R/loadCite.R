#' @title Load CITE-seq data
#' @description This function loads CITE-seq data from specified directories and samples.
#' @param directory A character vector of directories containing the CITE-seq data.
#' @param sample A character vector of sample names corresponding to the directories.
#' @param batch_barcode A character vector of batch barcodes to filter the data.
#' @return A SingleCellExperiment object containing the loaded CITE-seq data.
#' @export

loadCite <- function(directory, sample, batch_barcode = NULL) {
  sce_list <- list()
  for (i in 1:length(directory)) {
    dir <- directory[i]
    name <- sample[i]
    temp <- DropletUtils::read10xCounts(dir, sample.names = name, col.names = TRUE)
    temp$Barcode <- paste0(temp$Barcode, "-1")
    temp$batch_barcode <- paste0(temp$Sample, "_", temp$Barcode)
    colnames(temp) <- temp$batch_barcode

    if (!is.null(batch_barcode)) {
      temp <- temp[ ,colnames(temp) %in% batch_barcode]
    }

    temp <- temp[!rownames(temp) %in% "unmapped", ]
    assay(temp, "counts") <- as.matrix(assay(temp, "counts"))
    sce_list[[name]] <- temp
    message(paste0("Finish reading data: ", name))
  }

  if (length(sample) == 1){
    message("Number of samples is 1, returning SingleCellExperiment object...")
    sce <- sce_list[[1]]
  } else {
    message(paste0("Number of samples is ", length(sample), ", combining CITE-seq data..."))
    
    # Get common features across all samples (intersection)
    all_features <- lapply(sce_list, rownames)
    common_features <- Reduce(intersect, all_features)
    message(paste0("Found ", length(common_features), " common features across all samples"))
    
    # Subset each SCE to common features and combine count matrices
    combined_counts_list <- list()
    combined_coldata_list <- list()
    combined_rowdata_list <- list()
    
    for (i in 1:length(sce_list)) {
      sce_temp <- sce_list[[i]]
      # Subset to common features
      sce_temp <- sce_temp[common_features, ]
      # Extract counts matrix
      combined_counts_list[[i]] <- counts(sce_temp)
      # Extract and modify colData
      temp_coldata <- colData(sce_temp)
      temp_coldata$batch <- sample[i]
      combined_coldata_list[[i]] <- temp_coldata
      # Store rowData (should be the same for common features)
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
    
    message("Done combining CITE-seq data...")
  }
  rowData(sce)$name <- rownames(sce)
  rowdata <- strsplit(rowData(sce)$name, "-\\s*(?=[^-]+$)", perl=TRUE)
  rowdata <- as.data.frame(t(as.data.frame(rowdata)))
  colnames(rowdata) <- c("cite", "barcode")
  rowdata$cite <- paste0("cite_", rowdata$cite)
  rowData(sce) <- cbind(rowData(sce), rowdata)
  rownames(sce) <- rowData(sce)$cite
  sce$batch_barcode <- colnames(sce)
  return(sce)
}
