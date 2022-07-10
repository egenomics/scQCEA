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
