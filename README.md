# bioscrna

R package for a Bioconductor-based single-cell RNA-seq analysis pipeline: 10X loading, QC, normalization, dimension reduction, batch correction, clustering, marker detection, and plotting helpers for `SingleCellExperiment` objects.

## Citation

If you use bioscrna, please cite:

Wu, X., Zhang, M., Moon, J.S. *et al.* Senescent-like activated naive B cells contribute to progression to clinical rheumatoid arthritis. *Nat Aging* **6**, 2143–2163 (2026). https://doi.org/10.1038/s43587-026-01236-6

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

Copyright (C) 2026 Mengrui Zhang.

bioscrna is released under the [GNU General Public License v3.0](LICENSE) (GPL-3). You may use, modify, and distribute it, including for commercial purposes, if you keep the copyright notice and release distributed derivatives under GPL-3 as well. That requires credit and blocks closed-source redistribution of this code.
