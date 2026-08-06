# bioscrna

R package for a Bioconductor-based single-cell RNA-seq analysis pipeline: 10X loading, QC, normalization, dimension reduction, batch correction, clustering, marker detection, and plotting helpers for `SingleCellExperiment` objects.

## Citation

If you use bioscrna, please cite:

Xiaohao Wu, Mengrui Zhang, Jae-Seung Moon, Eun Kyung Song, Joshua C. Abrams, Laura S. van Dam, Orr Sharpe, Marie Feser, Laurie Moss, Peggy P. Ho, Melanie H. Smith, Laura T. Donlin, Jane H. Buckner, Eddie A. James, Gary S. Firestein, Yuko Okamoto, Tobiaz V. Lanz, Eric Meffre, V. Michael Holers, Kevin D. Deane, William H. Robinson. Senescent Activated Naive B Cells Promote Anti-Citrullinated Antigen T Cell Responses and the Transition to Clinical Rheumatoid Arthritis. bioRxiv. 2025. doi:[10.1101/2025.10.14.682430](https://doi.org/10.1101/2025.10.14.682430).

## Installation

```r
# install.packages("remotes")
remotes::install_github("zhangzmr/bioscrna")
```

For local development from a clone:

```r
# install.packages("devtools")
devtools::install("path/to/bioscrna")
```

## Optional dependency

Batch correction with Harmony requires the [harmony](https://github.com/immunogenomics/harmony) package:

```r
install.packages("harmony")
```

## Example workflow

```r
library(bioscrna)

# sce <- loading_sce(directory = dirs, sample = samples)
# plotQC(sce, file = "qc_before.pdf")
# sce <- qualityfilter(sce)
# sce <- run_rna(sce, batch = "batch", merge_method = c("FastMNN", "Harmony"))
# res <- clustering(sce, k_range = 5:20, dimred = "UMAP_corrected")
# sce <- res$sce
# plotCluster(sce, dimred = "UMAP_corrected", file = "clusters.pdf")
# markers <- clusterMarker(sce, top = 5, cluster = sce$k.10_cluster.fun.louvain)
```

## License

MIT © Mengrui Zhang
