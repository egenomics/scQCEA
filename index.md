## **scQCEA**

### Introduction 
This documentation gives an introduction and usage manual of scQCEA (A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data), for annotating and interpreting scRNA-Seq data.
<br />
scQCEA is an R package for annotation and quality control report of scRNA-Seq profiles. It generates an interactive report of quality control metrics which allows visual evaluation of QC metrics, objective selection of insightful optimal cluster numbers and discrimination between true variation and background noise.
<br />
The easiest way to generate an interactive summary QC report is to run the `RUN_ME.R` script from the RStudio. The required inputs are a gene-cell count matrix, feature-barcode matrices, and tSNE and UMAP projections from 10X CellRanger count.

### Easy Installation
1. Install the R (https://cran.r-project.org/)
<br />
2. Install the free version of rStudio (https://www.rstudio.com/products/rstudio/download/)
<br />
3. Download scQCEA from GitHub [LINK](https://github.com/isarnassiri/scQCEA/)
[Figure 1](https://github.com/isarnassiri/scQCEA/tree/gh-pages/Download_Github.png)
<br />
4. To install scQCEA, run the `RUN_ME.R` script from the RStudio. All dependency packages automatically will be downloaded, installed and loaded from CRAN-like repositories. The following versions of its dependencies are compatible with the library:

**Required version of packages in R:**

```markdown
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
```

**Environment:** 
<br />
We only tested scQCEA in the R version 4.1.1 (2021-08-10) environment. You need to have root permission for this distribution, including the installation of any package.

### Install from Source Code
Alternatively, you can download the source codes from CRAN and install libraries using the terminal as follows:

```markdown
# Required packages in R:

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
```

### Manual
It is easy to create an interactive QC report for those who possess little or no programming language skills. To run and generate an interactive QC report on your computer please open the `RUN_ME.R` file using rStudio, and click on the "Run" icon. An interactive QC report automatically will be generated in one HTML file. The scQCEA generates a QC report as an HTML file including four sections: experimental workflow, data processing workflow, sample information and QC metrics, data analysis and quality control.
<br />
By default, the HTML report will be written in /Outputs directory named `CLICK_ME.html`. You can open `CLICK_ME.html` without using rStudio/R. In addition, you can find a zip file in the /Outputs directory which is particularly useful to share or store the QC reports. 

### Input Data
As input, the scQCEA package expects raw count data from 10X CellRanger or other single-cell experiments and optional arguments such as appropriate organism. The gene-cell count data has the gene as a row (the gene name should be the human or mouse Ensembl gene ID) and the cell as a column. You can convert an HDF5 Feature-Barcode Matrix to a gene-cell count matrix using the cellranger mat2csv command provided by 10Xgenomics. The tSNE and UMAP projections are the outputs of dimensionality reduction analysis in CSV format.

### Cell Type Enrichment Analysis
Cell type annotation on scRNA-Seq data is a pre-step for generating an interactive QC report with scQCEA. This step requires some bioinformatics efforts, but there are a few good existing software to use.

**Recommended strategy for cell-type enrichment analysis:**
<br />
AUCell algorithm can be applied to score the activity of each reference gene set in each cell (Aibar, et al., 2017). It uses the area under the curve (AUC) to quantify the enrichment of an indicated reference gene set among the most highly expressed genes in each cell. The cell-type enrichment analysis only takes the count matrix from CellRanger count as the variable input.

### History
**Release v0.1.1 (04/07/2022)**
A completed version for all planned features.

### Quick Resources
Latest version on GitHub [LINK](https://github.com/isarnassiri/scQCEA/)

### Issue Reports
If you find any error or suspicious bug, we will appreciate your report. Please write them in the github issues: [LINK](https://github.com/isarnassiri/scQCEA/issues)

### References
Aibar, S., et al. SCENIC: single-cell regulatory network inference and clustering. Nature Methods 2017;14(11):1083-1086.
<br />
Fairfax, B.P., et al. Peripheral CD8+ T cell characteristics associated with durable responses to immune checkpoint blockade in patients with metastatic melanoma. Nature Medicine 2020;26(2):193-199.
<br />
Nassiri, I., Fairfax, B., Lee, A., Wu, Y., Buck, D., Piazza, P. scQCEA: A Framework for Annotation and Quality Control Report of Single-Cell RNA-Sequencing Data. 
