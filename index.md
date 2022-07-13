## **scQCEA**

### Introduction 
This documentation gives an introduction and usage manual of scQCEA (A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data), for annotating and interpreting scRNA-Seq data.
<br />
scQCEA is an R package for annotation and quality control report of scRNA-Seq profiles. It generates an interactive report of quality control metrics which allows visual evaluation of QC metrics, objective selection of insightful optimal cluster numbers and discrimination between true variation and background noise.
<br />
The easiest way to generate an interactive summary QC report is to run the `RUN_ME.R` script from the RStudio. The required inputs are a gene-cell count matrix, feature-barcode matrices, and tSNE and UMAP projections from 10X CellRanger count.

### Easy Installation
1. Install the R [(LINK)](https://cran.r-project.org/)
2. Install the free version of rStudio [(LINK)](https://www.rstudio.com/products/rstudio/download/)
3. Download scQCEA from GitHub [(LINK)](https://github.com/isarnassiri/scQCEA/), and unzip the folder

![Figure 1](/Download_Github.png)
4. To install scQCEA, run the `RUN_ME.R` script from the RStudio. All dependency packages automatically will be downloaded, installed and loaded from CRAN-like repositories. The following versions of its dependencies are compatible with the library:

**Required version of packages in R:**

```markdown
─ Session info ──────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.1.3 (2022-03-10)
 os       macOS Big Sur 11.6.3
 system   x86_64, darwin17.0
 ui       RStudio
 language (EN)
 collate  en_GB.UTF-8
 ctype    en_GB.UTF-8
 tz       Europe/London
 date     2022-07-13
 rstudio  2022.02.1+461 Prairie Trillium (desktop)
 pandoc   2.17.1.1 @ /Applications/RStudio.app/Contents/MacOS/quarto/bin/ (via rmarkdown)

─ Packages ──────────────────────────────────────────────────────────────────────────────────────
 package              * version  date (UTC) lib source
 annotate             * 1.72.0   2021-10-26 [1] Bioconductor
 AnnotationDbi        * 1.56.2   2021-11-09 [1] Bioconductor
 askpass                1.1      2019-01-13 [1] CRAN (R 4.1.0)
 assertthat             0.2.1    2019-03-21 [1] CRAN (R 4.1.0)
 AUCell               * 1.16.0   2021-10-26 [1] Bioconductor
 beachmat               2.10.0   2021-10-26 [1] Bioconductor
 Biobase              * 2.54.0   2021-10-26 [1] Bioconductor
 BiocGenerics         * 0.40.0   2021-10-26 [1] Bioconductor
 BiocManager          * 1.30.17  2022-04-22 [1] CRAN (R 4.1.2)
 BiocParallel           1.28.3   2021-12-09 [1] Bioconductor
 Biostrings             2.62.0   2021-10-26 [1] Bioconductor
 bit                    4.0.4    2020-08-04 [1] CRAN (R 4.1.0)
 bit64                  4.0.5    2020-08-30 [1] CRAN (R 4.1.0)
 bitops                 1.0-7    2021-04-24 [1] CRAN (R 4.1.0)
 blob                   1.2.3    2022-04-10 [1] CRAN (R 4.1.2)
 brio                   1.1.3    2021-11-30 [1] CRAN (R 4.1.0)
 bslib                  0.3.1    2021-10-06 [1] CRAN (R 4.1.0)
 bsselectR            * 0.1.0    2022-04-18 [1] Github (walkerke/bsselectR@c196e8f)
 cachem                 1.0.6    2021-08-19 [1] CRAN (R 4.1.0)
 callr                  3.7.0    2021-04-20 [1] CRAN (R 4.1.0)
 cli                    3.3.0    2022-04-25 [1] CRAN (R 4.1.2)
 cluster              * 2.1.2    2021-04-17 [1] CRAN (R 4.1.3)
 codetools              0.2-18   2020-11-04 [1] CRAN (R 4.1.3)
 colorspace             2.0-3    2022-02-21 [1] CRAN (R 4.1.2)
 crayon                 1.5.1    2022-03-26 [1] CRAN (R 4.1.2)
 crosstalk              1.2.0    2021-11-04 [1] CRAN (R 4.1.0)
 data.table           * 1.14.2   2021-09-27 [1] CRAN (R 4.1.0)
 DBI                    1.1.2    2021-12-20 [1] CRAN (R 4.1.0)
 DelayedArray           0.20.0   2021-10-26 [1] Bioconductor
 DelayedMatrixStats     1.16.0   2021-10-26 [1] Bioconductor
 desc                   1.4.1    2022-03-06 [1] CRAN (R 4.1.2)
 devtools             * 2.4.3    2021-11-30 [1] CRAN (R 4.1.0)
 DiagrammeR           * 1.0.9    2022-03-05 [1] CRAN (R 4.1.2)
 digest                 0.6.29   2021-12-01 [1] CRAN (R 4.1.0)
 doParallel             1.0.17   2022-02-07 [1] CRAN (R 4.1.2)
 downloadthis         * 0.3.1    2022-02-23 [1] CRAN (R 4.1.2)
 dplyr                * 1.0.9    2022-04-28 [1] CRAN (R 4.1.2)
 dqrng                  0.3.0    2021-05-01 [1] CRAN (R 4.1.0)
 DropletUtils         * 1.14.2   2022-01-09 [1] Bioconductor
 DT                   * 0.22     2022-03-28 [1] CRAN (R 4.1.2)
 edgeR                  3.36.0   2021-10-26 [1] Bioconductor
 ellipsis               0.3.2    2021-04-29 [1] CRAN (R 4.1.0)
 evaluate               0.15     2022-02-18 [1] CRAN (R 4.1.2)
 fansi                  1.0.3    2022-03-24 [1] CRAN (R 4.1.2)
 fastmap                1.1.0    2021-01-25 [1] CRAN (R 4.1.0)
 foreach                1.5.2    2022-02-02 [1] CRAN (R 4.1.2)
 fs                     1.5.2    2021-12-08 [1] CRAN (R 4.1.0)
 generics               0.1.2    2022-01-31 [1] CRAN (R 4.1.2)
 GenomeInfoDb         * 1.30.1   2022-01-30 [1] Bioconductor
 GenomeInfoDbData       1.2.7    2022-05-10 [1] Bioconductor
 GenomicRanges        * 1.46.1   2021-11-18 [1] Bioconductor
 GEOquery             * 2.62.2   2022-01-11 [1] Bioconductor
 ggplot2              * 3.3.6    2022-05-03 [1] CRAN (R 4.1.2)
 glue                   1.6.2    2022-02-24 [1] CRAN (R 4.1.2)
 graph                * 1.72.0   2021-10-26 [1] Bioconductor
 gridBase               0.4-7    2014-02-24 [1] CRAN (R 4.1.0)
 GSEABase             * 1.56.0   2021-10-26 [1] Bioconductor
 gtable                 0.3.0    2019-03-25 [1] CRAN (R 4.1.0)
 HDF5Array              1.22.1   2021-11-14 [1] Bioconductor
 hms                    1.1.1    2021-09-26 [1] CRAN (R 4.1.0)
 htmltools              0.5.2    2021-08-25 [1] CRAN (R 4.1.0)
 htmlwidgets            1.5.4    2021-09-08 [1] CRAN (R 4.1.0)
 httpuv                 1.6.5    2022-01-05 [1] CRAN (R 4.1.2)
 httr                   1.4.3    2022-05-04 [1] CRAN (R 4.1.2)
 IRanges              * 2.28.0   2021-10-26 [1] Bioconductor
 iterators              1.0.14   2022-02-05 [1] CRAN (R 4.1.2)
 jquerylib              0.1.4    2021-04-26 [1] CRAN (R 4.1.0)
 jsonlite               1.8.0    2022-02-22 [1] CRAN (R 4.1.2)
 kableExtra           * 1.3.4    2021-02-20 [1] CRAN (R 4.1.0)
 KEGGREST               1.34.0   2021-10-26 [1] Bioconductor
 knitr                  1.39     2022-04-26 [1] CRAN (R 4.1.2)
 later                  1.3.0    2021-08-18 [1] CRAN (R 4.1.0)
 lattice                0.20-45  2021-09-22 [1] CRAN (R 4.1.3)
 lazyeval               0.2.2    2019-03-15 [1] CRAN (R 4.1.0)
 lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.0)
 limma                  3.50.3   2022-04-07 [1] Bioconductor
 locfit                 1.5-9.5  2022-03-03 [1] CRAN (R 4.1.2)
 magrittr               2.0.3    2022-03-30 [1] CRAN (R 4.1.2)
 Matrix               * 1.4-0    2021-12-08 [1] CRAN (R 4.1.3)
 MatrixGenerics       * 1.6.0    2021-10-26 [1] Bioconductor
 matrixStats          * 0.62.0   2022-04-19 [1] CRAN (R 4.1.2)
 memoise                2.0.1    2021-11-26 [1] CRAN (R 4.1.0)
 mime                   0.12     2021-09-28 [1] CRAN (R 4.1.0)
 munsell                0.5.0    2018-06-12 [1] CRAN (R 4.1.0)
 NMF                  * 0.24.0   2022-03-29 [1] CRAN (R 4.1.2)
 pdftools             * 3.3.0    2022-07-07 [1] CRAN (R 4.1.2)
 pillar                 1.7.0    2022-02-01 [1] CRAN (R 4.1.2)
 pkgbuild               1.3.1    2021-12-20 [1] CRAN (R 4.1.0)
 pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.1.0)
 pkgload                1.2.4    2021-11-30 [1] CRAN (R 4.1.0)
 pkgmaker             * 0.32.2   2020-10-20 [1] CRAN (R 4.1.0)
 plotly               * 4.10.0   2021-10-09 [1] CRAN (R 4.1.0)
 plyr                   1.8.7    2022-03-24 [1] CRAN (R 4.1.2)
 png                  * 0.1-7    2013-12-03 [1] CRAN (R 4.1.0)
 prettyunits            1.1.1    2020-01-24 [1] CRAN (R 4.1.0)
 processx               3.5.3    2022-03-25 [1] CRAN (R 4.1.2)
 promises               1.2.0.1  2021-02-11 [1] CRAN (R 4.1.0)
 ps                     1.7.0    2022-04-23 [1] CRAN (R 4.1.2)
 purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.1.0)
 qpdf                   1.2.0    2022-05-29 [1] CRAN (R 4.1.2)
 R.methodsS3          * 1.8.1    2020-08-26 [1] CRAN (R 4.1.0)
 R.oo                 * 1.24.0   2020-08-26 [1] CRAN (R 4.1.0)
 R.utils              * 2.11.0   2021-09-26 [1] CRAN (R 4.1.0)
 R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.0)
 RColorBrewer           1.1-3    2022-04-03 [1] CRAN (R 4.1.2)
 Rcpp                   1.0.8.3  2022-03-17 [1] CRAN (R 4.1.2)
 RCurl                  1.98-1.6 2022-02-08 [1] CRAN (R 4.1.2)
 readr                * 2.1.2    2022-01-30 [1] CRAN (R 4.1.2)
 registry             * 0.5-1    2019-03-05 [1] CRAN (R 4.1.0)
 remotes                2.4.2    2021-11-30 [1] CRAN (R 4.1.0)
 reshape2               1.4.4    2020-04-09 [1] CRAN (R 4.1.0)
 rhdf5                  2.38.1   2022-03-10 [1] Bioconductor
 rhdf5filters           1.6.0    2021-10-26 [1] Bioconductor
 Rhdf5lib               1.16.0   2021-10-26 [1] Bioconductor
 rlang                  1.0.2    2022-03-04 [1] CRAN (R 4.1.2)
 rmarkdown            * 2.13     2022-03-10 [1] CRAN (R 4.1.3)
 rngtools             * 1.5.2    2021-09-20 [1] CRAN (R 4.1.0)
 rprojroot              2.0.3    2022-04-02 [1] CRAN (R 4.1.2)
 RSQLite                2.2.14   2022-05-07 [1] CRAN (R 4.1.2)
 rstudioapi           * 0.13     2020-11-12 [1] CRAN (R 4.1.0)
 rvest                  1.0.2    2021-10-16 [1] CRAN (R 4.1.0)
 S4Vectors            * 0.32.4   2022-03-29 [1] Bioconductor
 sass                   0.4.1    2022-03-23 [1] CRAN (R 4.1.2)
 scales                 1.2.0    2022-04-13 [1] CRAN (R 4.1.2)
 scuttle                1.4.0    2021-10-26 [1] Bioconductor
 sessioninfo            1.2.2    2021-12-06 [1] CRAN (R 4.1.0)
 shiny                  1.7.1    2021-10-02 [1] CRAN (R 4.1.0)
 SingleCellExperiment * 1.16.0   2021-10-26 [1] Bioconductor
 sparseMatrixStats      1.6.0    2021-10-26 [1] Bioconductor
 stringi                1.7.6    2021-11-29 [1] CRAN (R 4.1.0)
 stringr              * 1.4.0    2019-02-10 [1] CRAN (R 4.1.0)
 SummarizedExperiment * 1.24.0   2021-10-26 [1] Bioconductor
 svglite                2.1.0    2022-02-03 [1] CRAN (R 4.1.2)
 systemfonts            1.0.4    2022-02-11 [1] CRAN (R 4.1.2)
 testthat               3.1.4    2022-04-26 [1] CRAN (R 4.1.2)
 tibble                 3.1.7    2022-05-03 [1] CRAN (R 4.1.2)
 tidyr                  1.2.0    2022-02-01 [1] CRAN (R 4.1.2)
 tidyselect             1.1.2    2022-02-21 [1] CRAN (R 4.1.2)
 tinytex                0.38     2022-03-29 [1] CRAN (R 4.1.2)
 tzdb                   0.3.0    2022-03-28 [1] CRAN (R 4.1.2)
 usethis              * 2.1.5    2021-12-09 [1] CRAN (R 4.1.0)
 utf8                   1.2.2    2021-07-24 [1] CRAN (R 4.1.0)
 vctrs                  0.4.1    2022-04-13 [1] CRAN (R 4.1.2)
 viridisLite            0.4.0    2021-04-13 [1] CRAN (R 4.1.0)
 visNetwork             2.1.0    2021-09-29 [1] CRAN (R 4.1.0)
 webshot                0.5.3    2022-04-14 [1] CRAN (R 4.1.2)
 withr                  2.5.0    2022-03-03 [1] CRAN (R 4.1.2)
 xfun                   0.30     2022-03-02 [1] CRAN (R 4.1.2)
 XML                  * 3.99-0.9 2022-02-24 [1] CRAN (R 4.1.2)
 xml2                   1.3.3    2021-11-30 [1] CRAN (R 4.1.0)
 xtable                 1.8-4    2019-04-21 [1] CRAN (R 4.1.0)
 XVector                0.34.0   2021-10-26 [1] Bioconductor
 yaml                   2.3.5    2022-02-21 [1] CRAN (R 4.1.2)
 zip                  * 2.2.0    2021-05-31 [1] CRAN (R 4.1.0)
 zlibbioc               1.40.0   2021-10-26 [1] Bioconductor

 [1] /Library/Frameworks/R.framework/Versions/4.1/Resources/library

─────────────────────────────────────────────────────────────────────────────────────────────────
   
```

**Environment:** 
<br />
We only tested scQCEA in the R version 4.1.3 (2022-03-10) environment. You need to have root permission for this distribution, including the installation of any package.

### Install from Source Code
Alternatively, you can download the source codes from CRAN and install libraries using the terminal.

### Manual
It is easy to create an interactive QC report for those who possess little or no programming language skills. To run and generate an interactive QC report on your computer please open the `RUN_ME.R` file using rStudio, and click on the "Run" icon. An interactive QC report automatically will be generated in one HTML file, including four sections: experimental workflow, data processing workflow, sample information and QC metrics, data analysis and quality control (Fig. 2).


![Figure 2](/Figure_1.png)

Experimental workflow describes scRNA-seq transcriptome processing and sequencing platform. Data processing workflow presents an analysis pipeline to process data, including aligning reads, generating feature-barcode matrices, and other secondary analyses. Samples information and QC metrics provide tables of metadata and QC, listing a variety of metrics per application. Data analysis and quality control present projection of transcriptionally and functionally distinct clusters, highlighted by cell type group, including UMAP and t-SNE plots. Diagnostic plots provide technical features, including the distribution of non-duplicate reads with mapping quality per barcode.
<br />
By default, the HTML report will be written in /Outputs directory named `CLICK_ME.html`. You can open `CLICK_ME.html` without using rStudio/R. In addition, you can find a zip file in the /Outputs directory which is particularly useful to share or store the QC reports. The content of the "Data processing Workflow" section is automatically adjusted based on the type of application (s) and the "Library Type" column in "samples.metadata" file.
<br />

### Input Data
As input, the scQCEA package expects raw count data from 10X CellRanger (/outs/metrics_summary.csv) or other single-cell experiments and a metadata files arguments such as appropriate organism. The gene-cell count data has the gene as a row (the gene name should be the human or mouse Ensembl gene ID) and the cell as a column. You can convert an HDF5 Feature-Barcode Matrix to a gene-cell count matrix using the cellranger mat2csv command provided by 10Xgenomics.  The tSNE and UMAP projections are the outputs of dimensionality reduction analysis in CSV format.

### Cell Type Enrichment Analysis
Cell type annotation on scRNA-Seq data is a pre-step for generating an interactive QC report with scQCEA. This step requires some bioinformatics efforts, but there are a few good existing software to use.

**Recommended strategy for cell-type enrichment analysis:**
<br />
AUCell algorithm can be applied to score the activity of each reference gene set in each cell (Aibar, et al., 2017). It uses the area under the curve (AUC) to quantify the enrichment of an indicated reference gene set among the most highly expressed genes in each cell. The required inputs are a gene-cell count matrix, feature-barcode matrices, and tSNE and UMAP projections from 10X CellRanger count.

### History
**Release v0.1.1 (04/07/2022)**
A completed version for all planned features.

### Quick Resources
Latest version on GitHub [(LINK)](https://github.com/isarnassiri/scQCEA/)

### Issue Reports
If you find any error or suspicious bug, we will appreciate your report. Please write them in the github issues [(LINK)](https://github.com/isarnassiri/scQCEA/issues)

### References
1. Aibar, S., et al. SCENIC: single-cell regulatory network inference and clustering. Nature Methods 2017;14(11):1083-1086.
2. Fairfax, B.P., et al. Peripheral CD8+ T cell characteristics associated with durable responses to immune checkpoint blockade in patients with metastatic melanoma. Nature Medicine 2020;26(2):193-199.
3. Nassiri, I., Fairfax, B., Lee, A., Wu, Y., Buck, D., Piazza, P. scQCEA: A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data. 
