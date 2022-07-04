## scQCEA

### Introduction 

This documentation gives an introduction and usage manual of scQCEA (A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data), for annotating and interpreting scRNA-Seq data.

scQCEA is an R package for annotation and quality control report of scRNA-Seq profiles. It generates an interactive report of quality control metrics which allows visual evaluation of QC metrics, objective selection of insightful optimal cluster numbers and discrimination between true variation and background noise.

The easiest way to generate an interactive summary QC report is to run the RUN_ME.R script from the RStudio. The required inputs are a gene-cell count matrix, feature-barcode matrices, and tSNE and UMAP projections from 10X CellRanger count.

### Easy Installation

1. Install the R (https://cran.r-project.org/)

2. Install the free version of rStudio (https://www.rstudio.com/products/rstudio/download/)

3. To install scQCEA, run the RUN_ME.R script from the RStudio. All dependency packages automatically will be downloaded, installed and loaded from CRAN-like repositories. The following versions of its dependencies are compatible with the library:

**Required packages in R:**

readr_2.1.0       
dplyr_1.0.7        
data.table_1.14.2 
ggplot2_3.3.5      
downloadthis_0.2.1
DT_0.19            
bsselectR_0.1.0    
stringr_1.4.0     
rstudioapi_0.13    
rmarkdown_2.11  
kableExtra_1.3.4
R.utils_2.11.0 

**Environment:** 

We only tested scQCEA in the R version 4.1.1 (2021-08-10) environment. You need to have root permission for this distribution, including the installation of any package.

### Install from source code

Alternatively, you can download the source codes from CRAN and install libraries using the terminal as follows:

library(stringr)
library(bsselectR)
library(kableExtra)
library(DT)
library(downloadthis)
library(ggplot2)
library(data.table)
library(dplyr)
library(readr)
library(rmarkdown)
library(R.utils)
library(rstudioapi)


### Manual

It is easy to create an interactive QC report for those who possess little or no programming language skills. To run and generate an interactive QC report on your computer please open the "RUN_ME.R" file using rStudio, and click on the "Run" icon. An interactive QC report automatically will be generated in one HTML file. The scQCEA generates a QC report as an HTML file including four sections: experimental workflow, data processing workflow, samples information and QC metrics, data analysis and quality control.

By default, the HTML report will be written in /Outputs directory named CLICK_ME.html. You can open CLICK_ME.html without using rStudio/R. In addition, you can find a zip file in the /Outputs directory which is particularly useful to share or store the QC reports. 


### Input data

As input, the scQCEA package expects raw count data from 10X CellRanger or other single cell experiments, and optional arguments such as appropriate organism. The gene-cell count data has the gene as row (the gene name should be the human or mouse Ensembl gene ID) and the cell as column. You can convert a HDF5 Feature-Barcode Matrix to gene-cell count matrix using the cellranger mat2csv command provided by 10Xgenomics. The tSNE and UMAP projections are the outputs of dimensionality reduction analysis in csv format.

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/isarnassiri/scQCEA/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.

### Quick Resources
Latest version on GitHub https://github.com/isarnassiri/scQCEA/

Scripts for simulation https://github.com/single-cell-genetics/vireo/tree/master/simulate

All releases https://pypi.org/project/vireoSNP/#history

### Issue reports
If you find any error or suspicious bug, we will appreciate your report. Please write them in the github issues:(https://github.com/isarnassiri/scQCEA/issues)

### References
scQCEA: A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data. Isar Nassiri, Benjamin Fairfax, Angela Lee, Yanxia Wu, David Buck, Paolo Piazza.


