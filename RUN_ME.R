#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

# Install and load R package
if(length(list.of.packages[!("rstudioapi" %in% installed.packages()[,"Package"])])) install.packages(list.of.packages[!("rstudioapi" %in% installed.packages()[,"Package"])], repos = "http://cran.us.r-project.org")
library("rstudioapi")
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/Scripts/')); 

source("Generate_Interactive_QC_Report.R")

# Generate an "Interactive QC Report"
Generate_Interactive_QC_Report()

############################################################ 
# 'Find the "Interactive QC Report" in the Outputs/ folder #
############################################################
