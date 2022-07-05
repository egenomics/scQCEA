scQCEA
==========
* [Introduction](#introduction)
* [Installation](#Installation)
* [Usage](#Usage)
* [Citation](#citation)
<a name="introduction"/>

### Introduction

scQCEA is an R package for annotation and quality control report of scRNA-Seq profiles, which performs a probabilistic assignment of the reference cell types to identify clusters, before downstream analysis such as gene network inference. scQCEA provides automated cell type annotation on scRNA-seq data and identifies differential patterns in gene expression. scQCEA generates an interactive report of quality control metrics which allows visual evaluation of QC metrics, objective selection of insightful optimal cluster numbers and discrimination between true variation and background noise. 

Please see the [manual](https://isarnassiri.github.io/scQCEA/) for the usage of scRNABatchQC and the explanation of the HTML report.

<a name="installation"/>

### Installation
1. Install the R [(LINK)](https://cran.r-project.org/)
2. Install the free version of rStudio [(LINK)](https://www.rstudio.com/products/rstudio/download/)
3. Download scQCEA from GitHub [(LINK)](https://github.com/isarnassiri/scQCEA/), and unzip the folder
4. To install scQCEA, run the `RUN_ME.R` script from the RStudio. All dependency packages automatically will be downloaded, installed and loaded from CRAN-like repositories.

```{r,eval=FALSE}

#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

# Install and load R package
new.packages <- list.of.packages[!("rstudioapi" %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
library("rstudioapi") 
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/Scripts/')); 

source("Generate_Interactive_QC_Report.R")

# Generate an "Interactive QC Report"
Generate_Interactive_QC_Report()

############################################################ 
# 'Find the "Interactive QC Report" in the Outputs/ folder #
############################################################
```

### Usage

It is easy to create an interactive QC report for those who possess little or no programming language skills. To run and generate an interactive QC report on your computer please open the `RUN_ME.R` file using rStudio, and click on the "Run" icon. An interactive QC report automatically will be generated in one HTML file, including four sections: experimental workflow, data processing workflow, sample information and QC metrics, data analysis and quality control. 

By default, the HTML report will be written in /Outputs directory named `CLICK_ME.html`. You can open `CLICK_ME.html` without using rStudio/R. In addition, you can find a zip file in the /Outputs directory which is particularly useful to share or store the QC reports. The content of the "Data processing Workflow" section is automatically adjusted based on the type of application (s) and the "Library Type" column in "samples.metadata" file.

### Citation

Isar Nassiri, Benjamin Fairfax, Angela Lee, Yanxia Wu, David Buck, Paolo Piazza. scQCEA: A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data. 
