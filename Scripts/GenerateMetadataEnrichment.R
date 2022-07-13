
GenerateMetadataEnrichment <- function(){
  
# Install and load R packages
list.of.packages <- c("data.table", "rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
library(data.table); library(rstudioapi);

# getActiveDocumentContext shows the location of RNN_ME
setwd("~/"); setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/Inputs/'));
SamplesMetadata = fread('samples.metadata', stringsAsFactors = F, header = T)

# ---------------------------------- example input
project_id <- unique(SamplesMetadata$`Project Number`)
backend.data.dir <- "Inputs/ReferenceGeneSets/"

if(unique(SamplesMetadata$Genome) == 'GRCh38'){ organism <- "hsapiens" }
if(unique(SamplesMetadata$Genome) == 'GRCm38'){ organism <- "mmusculus" }
if(unique(SamplesMetadata$Genome) == 'GRCh38-premrna'){ organism <- "hsapiens" }
if(unique(SamplesMetadata$Genome) == 'GRCm38-premrna'){ organism <- "mmusculus" }

list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}

#list files recursive up to a certain level in R
dirs <- list.dirs.depth.n(".", n = 2)
#Include pattern in list.dirs
dirs <- grep("10X-gex", dirs, value = TRUE)

if(length(grep('./10X-gex-grouped/', dirs))>0)
{
  input.dir = dirs[grep('./10X-gex-grouped/', dirs)]
  input.dir = gsub('/./', '/', input.dir)
  
  OP = if(organism == "hsapiens"){paste0(backend.data.dir, 'human')}else{paste0(backend.data.dir, 'mouse')}
  
  guide_file = data.frame(
    project_id = replicate(length(input.dir),project_id),
    sample_id = gsub('.*\\/','',input.dir),
    inpute_dir = input.dir,
    backend_data_dir = replicate(length(input.dir),OP),
    organism = replicate(length(input.dir),organism))
  
  write.table(guide_file, 'gex_grouped_aggregation', quote = F, col.names = F, row.names = F, sep = '\t')
}

if(length(grep('./10X-gex/', dirs))>0)
{
  input.dir = dirs[grep('./10X-gex/', dirs)]
  input.dir = gsub('/./', '/', input.dir)
  
  OP = if(organism == "hsapiens"){paste0(backend.data.dir, 'human')}else{paste0(backend.data.dir, 'mouse')}
  
  guide_file = data.frame(
    project_id = replicate(length(input.dir),project_id),
    sample_id = gsub('.*\\/','',input.dir),
    inpute_dir = input.dir,
    backend_data_dir = replicate(length(input.dir),OP),
    organism = replicate(length(input.dir),organism))
  
  write.table(guide_file, 'gex_ungrouped_aggregation', quote = F, col.names = F, row.names = F, sep = '\t')
}
}
