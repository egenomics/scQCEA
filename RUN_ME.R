#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

##### Install and load R packages #####
if(length(c("devtools", "rstudioapi")[!(c("devtools", "rstudioapi") %in% installed.packages()[,"Package"])])) install.packages(c("devtools", "rstudioapi")[!(c("devtools", "rstudioapi") %in% installed.packages()[,"Package"])], repos = "http://cran.us.r-project.org", update = FALSE, dependencies = T); library("devtools", quietly = TRUE); library("rstudioapi", quietly = TRUE);

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
