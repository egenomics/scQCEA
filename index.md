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
 AUCell               * 1.16.0   2021-10-26 [1] Bioconductor
 Biobase              * 2.54.0   2021-10-26 [1] Bioconductor
 BiocGenerics         * 0.40.0   2021-10-26 [1] Bioconductor
 BiocManager          * 1.30.17  2022-04-22 [1] CRAN (R 4.1.2)
 bsselectR            * 0.1.0    2022-04-18 [1] Github (walkerke/bsselectR@c196e8f)
 cluster              * 2.1.2    2021-04-17 [1] CRAN (R 4.1.3)
 data.table           * 1.14.2   2021-09-27 [1] CRAN (R 4.1.0)
 devtools             * 2.4.3    2021-11-30 [1] CRAN (R 4.1.0)
 DiagrammeR           * 1.0.9    2022-03-05 [1] CRAN (R 4.1.2)
 downloadthis         * 0.3.1    2022-02-23 [1] CRAN (R 4.1.2)
 dplyr                * 1.0.9    2022-04-28 [1] CRAN (R 4.1.2)
 DropletUtils         * 1.14.2   2022-01-09 [1] Bioconductor
 DT                   * 0.22     2022-03-28 [1] CRAN (R 4.1.2)
 GenomeInfoDb         * 1.30.1   2022-01-30 [1] Bioconductor
 GenomicRanges        * 1.46.1   2021-11-18 [1] Bioconductor
 GEOquery             * 2.62.2   2022-01-11 [1] Bioconductor
 ggplot2              * 3.3.6    2022-05-03 [1] CRAN (R 4.1.2)
 graph                * 1.72.0   2021-10-26 [1] Bioconductor
 GSEABase             * 1.56.0   2021-10-26 [1] Bioconductor
 IRanges              * 2.28.0   2021-10-26 [1] Bioconductor
 kableExtra           * 1.3.4    2021-02-20 [1] CRAN (R 4.1.0)
 Matrix               * 1.4-0    2021-12-08 [1] CRAN (R 4.1.3)
 MatrixGenerics       * 1.6.0    2021-10-26 [1] Bioconductor
 matrixStats          * 0.62.0   2022-04-19 [1] CRAN (R 4.1.2)
 NMF                  * 0.24.0   2022-03-29 [1] CRAN (R 4.1.2)
 pdftools             * 3.3.0    2022-07-07 [1] CRAN (R 4.1.2)
 pkgmaker             * 0.32.2   2020-10-20 [1] CRAN (R 4.1.0)
 plotly               * 4.10.0   2021-10-09 [1] CRAN (R 4.1.0)
 png                  * 0.1-7    2013-12-03 [1] CRAN (R 4.1.0)
 R.methodsS3          * 1.8.1    2020-08-26 [1] CRAN (R 4.1.0)
 R.oo                 * 1.24.0   2020-08-26 [1] CRAN (R 4.1.0)
 R.utils              * 2.11.0   2021-09-26 [1] CRAN (R 4.1.0)
 readr                * 2.1.2    2022-01-30 [1] CRAN (R 4.1.2)
 registry             * 0.5-1    2019-03-05 [1] CRAN (R 4.1.0)
 rmarkdown            * 2.13     2022-03-10 [1] CRAN (R 4.1.3)
 rngtools             * 1.5.2    2021-09-20 [1] CRAN (R 4.1.0)
 rstudioapi           * 0.13     2020-11-12 [1] CRAN (R 4.1.0)
 S4Vectors            * 0.32.4   2022-03-29 [1] Bioconductor
 SingleCellExperiment * 1.16.0   2021-10-26 [1] Bioconductor
 stringr              * 1.4.0    2019-02-10 [1] CRAN (R 4.1.0)
 SummarizedExperiment * 1.24.0   2021-10-26 [1] Bioconductor
 usethis              * 2.1.5    2021-12-09 [1] CRAN (R 4.1.0)
 XML                  * 3.99-0.9 2022-02-24 [1] CRAN (R 4.1.2)
 zip                  * 2.2.0    2021-05-31 [1] CRAN (R 4.1.0)

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
