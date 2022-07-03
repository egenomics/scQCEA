
Generate_Interactive_QC_Report <- function(){
  
  # Install and load R packages
  list.of.packages <- c("rmarkdown", "rstudioapi", "zip")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
  library(rmarkdown); library(rstudioapi); library(zip);
  
  # Create an "Interactive QC Report"
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); 
  render("SourceCode.Rmd", quiet = T);
  
  invisible(file.rename('SourceCode.html', 'CLICK_ME.html'));
  invisible(file.copy(from = 'CLICK_ME.html', to = "Outputs/CLICK_ME.html"));
  invisible(file.remove(paste0(getwd(), '/CLICK_ME.html')));
  PInf = fread('Inputs/PInf.txt', stringsAsFactors = F, header = F)
  getwd()
  setwd('Outputs/'); 
  zip(zipfile = paste0('OGC_Interactive_QC_Report_', gsub('.*=','',PInf$V1[1]),'.zip'), files = c('CLICK_ME.html', 'Inputs'), recurse = TRUE, include_directories = TRUE);
  
}