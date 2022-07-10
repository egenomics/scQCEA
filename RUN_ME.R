#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

# Install and load R package
if(!("rstudioapi" %in% installed.packages()[,"Package"])) install.packages("rstudioapi", repos = "http://cran.us.r-project.org"); library("rstudioapi")
setwd("~/"); setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/Scripts/')); 

source("Generate_Interactive_QC_Report.R")

# Generate an "Interactive QC Report"
Generate_Interactive_QC_Report()

############################################################ 
#  Find the "Interactive QC Report" in the Outputs/ folder #
############################################################
