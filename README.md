scQCEA
==========
* [Introduction](#introduction)
* [Input data](#Input)
* [Cell type enrichment analysis](#CellTypeEnrichmentAnalysis)
* [Interactive QC report](#installation)
* [Citation](#citation)
<a name="introduction"/>

### Introduction

scQCEA is an R package for annotation and quality control report of scRNA-Seq profiles, which performs a probabilistic assignment of the reference cell types to identify clusters, before downstream analysis such as gene network inference. scQCEA provides automated cell type annotation on scRNA-seq data and identifies differential patterns in gene expression. scQCEA generates an interactive report of quality control metrics which allows visual evaluation of QC metrics, objective selection of insightful optimal cluster numbers and discrimination between true variation and background noise. Please see the [manual](https://isarnassiri.github.io/scQCEA/) for the usage of scRNABatchQC and the explanation of the HTML report.

<a name="installation"/>

### Installation
1. Install the R [(LINK)](https://cran.r-project.org/)
2. Install the free version of rStudio [(LINK)](https://www.rstudio.com/products/rstudio/download/)
3. Download scQCEA from GitHub [(LINK)](https://github.com/isarnassiri/scQCEA/), and unzip the folder
4. To install scQCEA, run the `RUN_ME.R` script from the RStudio. All dependency packages automatically will be downloaded, installed and loaded from CRAN-like repositories.

### Cell Type Enrichment analysis
Example of running CellTypeEnrichmentAnalysis function:

```{r,eval=FALSE}

library(scQCEA)
project_id='P180121'
gex_library_id='481207_03'
input.dir='/10X-gex/481207_03'
backend.data.dir='/Back_End_Data'
organism='hsapiens'
parent.dir='~/scQCEA/inst/extdata'
CellTypeEnrichmentAnalysis(project_id, gex_library_id, input.dir, backend.data.dir, organism, parent.dir)
```

### HTML report: interactive QC report

It is easy to create an interactive QC report for those who possess little or no programming language skills. To run and generate an interative QC report on your computer please:

1. Install the R (https://cran.r-project.org/)

2. Install the free version of rStudio (https://www.rstudio.com/products/rstudio/download/)

3. Open the "CLICK_ME.Rmd" file inside the folder [scQCEA_Interactive_Report] using rStudio, and click on the "Knite" blue icon. The user interface of the database will come up as a new window, and create an interactive QC report in one HTML file. By default, the HTML report will be written in your working directory named CLICK_ME.html. You can open and share CLICK_ME.html without using rStudio/R.

The scQCEA generates a QC report as an HTML file including four sections: experimental workflow, data processing workflow, samples information and QC metrics, data analysis and quality control.

### Citation

Isar Nassiri, Benjamin Fairfax, Angela Lee, Yanxia Wu, David Buck, Paolo Piazza. scQCEA: A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data. 
