
GenerateInteractiveQCReport <- function(){
  
  # Install and load R packages
  list.of.packages <- c("devtools", "rmarkdown", "rstudioapi", "zip")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org", update = FALSE, dependencies = T)

  suppressPackageStartupMessages({
    library(devtools, quietly = TRUE); library(rmarkdown, quietly = TRUE); library(rstudioapi, quietly = TRUE); library(zip, quietly = TRUE);
  })
  
  # Create an "Interactive QC Report"
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); 
  render("SourceCode.Rmd", quiet = T);
  
  invisible(file.rename('SourceCode.html', 'CLICK_ME.html'));
  invisible(file.copy(from = 'CLICK_ME.html', to = "Outputs/CLICK_ME.html"));
  invisible(file.remove(paste0(getwd(), '/CLICK_ME.html')));

  #--- Delete extra files
  setwd('Outputs/');

  list.dirs.depth.n <- function(p, n) {
    res <- list.dirs(p, recursive = FALSE)
    if (n > 1) {
      add <- list.dirs.depth.n(res, n-1)
      c(res, add)
    } else {
      res
    }
  }
  
  listdir = list.dirs.depth.n(".", n = 5)
  toMatch <- c('analysis', 'raw_feature_bc_matrix')
  matches <- unique(grep(paste(toMatch,collapse="|"), listdir, value=TRUE))
  invisible(unlink(matches, recursive = T, force = T));
  
  listfiles = list.files(pattern = 'read_count.csv', recursive = T)
  invisible(unlink(listfiles, recursive = T, force = T));
  
  listfiles = list.files(pattern = '*.pdf', recursive = T)
  invisible(unlink(listfiles, recursive = T, force = T));
  
  #--- generate zip file
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
  PInf = fread('Inputs/PInf.txt', stringsAsFactors = F, header = F);
  
  setwd('Outputs/');
  zip(zipfile = paste0('OGC_Interactive_QC_Report_', gsub('.*=','',PInf$V1[1]),'.zip'), files = c('CLICK_ME.html', 'Inputs'), recurse = TRUE, include_directories = TRUE);
  
}