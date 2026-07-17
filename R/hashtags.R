#' @title Hashtag demultiplexing
#' @description Demultiplexing of hashtag oligo (HTO) data using the hashedDrops function from the DropletUtils package.
#' @param rna_dirt Directory containing the RNA data.
#' @param cite_dirt Directory containing the CITE-seq data.
#' @param sample Sample name.
#' @param hashtag Vector of HTO names used in the experiment.
#' @param remove_doublet Logical, whether to remove doublet cells from the SingleCellExperiment object. Default is TRUE.
#' @return A SingleCellExperiment object with demultiplexed data and additional metadata columns.
#' @export
hashtags <- function(rna_dirt, cite_dirt, sample, hashtag, remove_doublet = TRUE) {
  sce <- loading_sce(directory = rna_dirt, sample = sample)
  cite <- loading_sce(directory = cite_dirt, sample = sample)
  batch_barcode <- intersect(sce$batch_barcode, cite$batch_barcode)
  sce <- sce[ ,sce$batch_barcode %in% batch_barcode]
  cite <- cite[ ,sce$batch_barcode]

  hto_counts <- assay(cite, "counts")
  hto_counts <- hto_counts[hashtag, ]
  demux <- data.frame(DropletUtils::hashedDrops(hto_counts))
  hashtag_mapping <- hashtag
  names(hashtag_mapping) <- seq_along(hashtag)
  demux$hashtag <- hashtag_mapping[as.character(demux$Best)]
  demux$hashlabels <- names(hashtag)[demux$Best]
  mucl <- demux[ ,c("Doublet", "Confident", "hashtag", "hashlabels")]
  colData(sce) <- cbind(colData(sce), mucl)
  message("Number of cells before removing doublets: ", ncol(sce))

  if (isTRUE(remove_doublet) && sum(!is.na(sce$Doublet)) > 0) {
    message("Filtering out doublet cells...")
    sce <- sce[ , sce$Doublet == FALSE]
    message("Number of cells after removing doublets: ", ncol(sce))
  } else {
    message("Doublet removal not performed.")
  }
  return(sce)
}





