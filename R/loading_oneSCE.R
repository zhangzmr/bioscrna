#' Load a single 10X Genomics dataset into a SingleCellExperiment object
#' This function reads the 10X Genomics data files (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
#' from the specified directory and creates a SingleCellExperiment object.
#' It also adds sample information to the colData and gene information to the rowData.
#' @param directory A character string specifying the directory containing the 10X Genomics data files.
#' @param sample A character string specifying the sample name to be added to the colData.
#' @return A SingleCellExperiment object containing the loaded data.
#' @importFrom Matrix readMM
#' @importFrom methods as
#' @importFrom utils read.delim
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
loading_oneSCE <- function(directory, sample) {
  # Read the 10X files manually
  message("Reading 10X matrix files...\n")
  counts <- readMM(file.path(directory, "matrix.mtx.gz"))
  features <- read.delim(file.path(directory, "features.tsv.gz"), header = FALSE, sep = "\t")
  barcodes <- read.delim(file.path(directory, "barcodes.tsv.gz"), header = FALSE, sep = "\t")
  counts <- as(counts, "CsparseMatrix")

  # Use gene symbols (column 2 if available, otherwise column 1) and make them unique to avoid duplicates
  if (ncol(features) >= 2) {
    gene_names <- features$V2
  } else {
    gene_names <- features$V1
  }
  gene_names <- make.unique(gene_names)

  # Set up matrix dimensions and names
  rownames(counts) <- gene_names
  colnames(counts) <- barcodes$V1

  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(assays = list(counts = counts))

  # Add sample information to colData
  sce$Sample <- sample
  sce$Barcode <- colnames(sce)
  sce$batch_barcode <- paste0(sce$Sample, "_", sce$Barcode)
  colnames(sce) <- sce$batch_barcode

  # Add gene information to rowData
  rowData(sce)$gene_id <- features$V1
  if (ncol(features) >= 2) {
    rowData(sce)$gene_symbol <- features$V2
  }
  if (ncol(features) >= 3) {
    rowData(sce)$gene_type <- features$V3
  }

  message("SingleCellExperiment object created successfully!\n")
  message("Dimensions:", nrow(sce), " genes x ", ncol(sce), " cells\n")
  return(sce)
}


