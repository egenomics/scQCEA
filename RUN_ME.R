#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

##### Install and load R packages #####
if(!("rstudioapi" %in% installed.packages()[,"Package"])) install.packages("rstudioapi", repos = "http://cran.us.r-project.org"); library("rstudioapi")

##### Generate an "Interactive QC Report" #####
setwd("~/"); setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/Scripts/')); 
source("GenerateInteractiveQCReport.R")
GenerateInteractiveQCReport()

##### Cell Type Enrichment Analysis #####
setwd("~/"); setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/Scripts/')); 
source("CellTypeEnrichment.R")
CellTypeEnrichment()

############################################################ 
#  Find the "Interactive QC Report" in the Outputs/ folder #
############################################################
