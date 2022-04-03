
#############################################################################################
###################################### Cell Type Enrichment #################################
#############################################################################################

#'@import AUCell
#'@import GSEABase
#'@import data.table
#'@import DT
#'@import NMF
#'@import plotly
#'@import GEOquery 
#'@import data.table ## 1.9.2 or 1.9.3
#'@import stringr
#'@import dplyr
#'@import argparser
#'@import optparse
#'@import dplyr
#'@import Matrix
#'@import DropletUtils
#'
#'@export
#'@name CellTypeEnrichmentAnalysis
#'@title Cell-Type Enrichment Analysis of scRNA-seq data
#'@description A Cell Type Annotation function for Single-Cell RNA-seq Data.
#'@author {Isar Nassiri}
#'@param project_id
#'Name of a folder which contains all input files
#'@param gex_library_id
#'Name of a folder which contains gene expression matrix
#'@param gex_library_id
#'Name of a folder which contains folders of gene expression matrix for several samples
#'@param backend.data.dir
#'Name of a folder which contains reference gene sets for different cell-types
#'@param organism
#'hsapiens or mmusculus
#'@return You can find the results as pdf and text files including cell type assignment matrix, tSNE, UMAP, and heatmap plots.
#'@examples
#'project_id='P180121'
#'gex_library_id=	'481207_03'
#'input.dir='/10X-gex/481207_03'
#'backend.data.dir='/Back_End_Data'
#'organism='hsapiens'
#'parent.dir='/users/buck/ani054/scQCEA/inst/extdata'
#'CellTypeEnrichmentAnalysis(project_id, gex_library_id, input.dir, backend.data.dir, organism, parent.dir)
#'@export

CellTypeEnrichmentAnalysis <- NULL
CellTypeEnrichmentAnalysis <- function(project_id, gex_library_id, input.dir, backend.data.dir, organism, parent.dir)
{
  
  # ----------------------------------
  project.id <- project_id
  gex.library.id <- gex_library_id
  input.dir <- input.dir
  parent.dir <- parent.dir
  backend.data.dir <- backend.data.dir
  organism <- organism
  
  setwd(parent.dir)
  
  if(length(grep('grouped', input.dir))==1) # catch integer(0)
  {output_dir = paste0(project_id, "/10X-gex-grouped/")
  }else{output_dir = paste0(project_id, "/10X-gex/")}
  
  # ---------------------------------- Validating arguments
  if (is.null(project_id)) {
    stop(paste0("Argument not provided. --project_id"))
  }
  if (!is.null(gex_library_id) && length(grep("[^a-zA-Z0-9_-]", gex_library_id)) > 0) {
    stop(paste0("Invalid characters detected. This option can only contain letters, numbers, dash (-) and underscore (_). --gex_library_id ", gex_library_id))
  }
  if (is.null(backend.data.dir)) {
    stop(paste0("Argument not provided or directory provided does not exist. --backend.data.dir ", backend.data.dir))
  }
  
  # ---------------------------------- Create output directory
  output.dir_perSample <- paste(output_dir, gex.library.id, sep="/")
  dir.create(output.dir_perSample, showWarnings=FALSE, recursive=TRUE) # recursive=TRUE will create parent directory first, if it doesn't exist.
  
  # ---------------------------------- read the scRNAseq profile
  TP_profile = fread(paste0(project.id, '/', input.dir, '/outs/read_count.csv'), stringsAsFactors = F)
  TP_profile = as.data.frame(TP_profile)
  row.names(TP_profile) = TP_profile$V1
  TP_profile_sub = TP_profile[,-1]
  colnames(TP_profile_sub) = gsub('-.*','',colnames(TP_profile_sub))
  
  # ---------------------------------- read repository of gene symbols and Ensemble transcript IDs
  if(organism == "hsapiens"){
    GTF = fread('ensembl_human.txt', stringsAsFactors = F, header = T)
    GTF = as.data.frame(GTF)
  }
  
  if(organism == "mmusculus"){
    GTF = fread('ensembl_moueSymID_humanSym.txt', stringsAsFactors = F, header = T)
    GTF = as.data.frame(GTF)
    GTF = GTF[,which(colnames(GTF) %in% c('gene_name', 'gene_id'))]
  }
  
  # ---------------------------------- add gene symbols to the gene expression profile and remove the gene IDs
  TP_profile_sub = data.frame(gene_id = row.names(TP_profile_sub), TP_profile_sub)
  TP_profile_sub = merge(GTF, TP_profile_sub, by = 'gene_id')
  
  TP_profile_sub$gene_name = make.names(TP_profile_sub$gene_name,unique=T)
  row.names(TP_profile_sub) = TP_profile_sub$gene_name
  TP_profile_sub = TP_profile_sub[,-c(1,2)]
  
  # ---------------------------------- enrichment input
  list = list.files('Back_End_Data', pattern = '.tsv')
  
  my_read_csv <- function(x) {
    out <- fread(x, stringsAsFactors = F, select = 'Gene')
    site <- gsub('blood_cell_category_rna_|blood_cell_category_rna_|.tsv.gz|_lineage|_Lineage|_Group','', x)
    site <- gsub('blood_cell_category_rna_|_Cell','', site)
    cbind(cellType=site, out)
  }
  
  repository <- lapply(paste('Back_End_Data',list,sep = '/'), my_read_csv)
  repository = do.call(rbind.data.frame, repository)
  
  # ---------------------------------- create a list of gene sets
  celltype = unique(repository$cellType)
  geneSets <- list()
  for(i in 1:length(celltype))
  {
    temp = repository[which(repository$cellType == celltype[i]),]
    geneSets[[i]] <- temp$Gene
  }
  
  names(geneSets) = gsub('.*/|.tsv', '', celltype)
  
  # ---------------------------------- AUCell
  set.seed(123)
  
  exprMatrix <- as.matrix(TP_profile_sub)
  cells_rankings <- AUCell_buildRankings(exprMatrix)
  dim(exprMatrix)
  # ---------------------------------- reports Genes in the gene sets NOT available in the dataset
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, nCores=1) 
  
  # ---------------------------------- cells assignment
  
  selectedThresholds <- NULL
  
  set.seed(123)
  par(mfrow=c(3,3)) # PROVIDE ENOUGH SPACE IN PLOTS ENVIRONMENT OF THE RSTUDIO
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE) 
  
  ## ------ check warnings
  # warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
  # warningMsg[which(warningMsg!="")]
  
  ## ------ export cell Assignment as text
  cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
  assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
  colnames(assignmentTable)[2] <- "geneSet"
  
  # ---------------------------------- save cell assignments as a heatmap
  assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
  
  set.seed(123)
  aheatmap(assignmentMat, scale="none", color="-RdYlBu2:100", legend=FALSE, filename = paste0(output.dir_perSample, '/', project.id, '_', gex.library.id, '_', 'Celltype_assignment_HeatMap.pdf'), height = 10, width = 8)
  
  # ---------------------------------- check if any cells assigned in more than one cell type
  # dim(assignmentTable)[1] > dim(exprMatrix)[2] # if it's true means if have duplications in enrichment results
  
  Freq_geneSet <- data.frame(table(assignmentTable[,"geneSet"]), stringsAsFactors = F)
  Freq_geneSet$Var1 = as.character(Freq_geneSet$Var1)
  colnames(Freq_geneSet)[1] = 'geneSet'
  
  # ---------------------------------- [Keep one assignment per cell - keep most frequently enriched one] Order data frame rows according to vector with specific order
  
  assignmentTable_rearranged = left_join(data.frame(geneSet=Freq_geneSet$geneSet),assignmentTable,by="geneSet")
  assignmentTable_dedup = assignmentTable_rearranged[!duplicated(assignmentTable_rearranged$cell),]
  
  # ---------------------------------- add number per enrichment cluster
  dim(assignmentTable_dedup)[1] > dim(exprMatrix)[2] # if it's TRUE means if have duplications in enrichment results
  
  Freq_assignmentTable_dedup = data.frame(table(assignmentTable_dedup[,"geneSet"]), stringsAsFactors = F)
  Freq_assignmentTable_dedup = Freq_assignmentTable_dedup[order(Freq_assignmentTable_dedup$Freq, decreasing = T),]
  Freq_assignmentTable_dedup$cluster = 1:dim(Freq_assignmentTable_dedup)[1]
  colnames(Freq_assignmentTable_dedup)[1] = 'geneSet'
  
  assignmentTable_dedup_ClusterNumber = left_join(assignmentTable_dedup, Freq_assignmentTable_dedup, by = 'geneSet')
  assignmentTable_dedup_ClusterNumber = assignmentTable_dedup_ClusterNumber[order(assignmentTable_dedup_ClusterNumber$Freq,decreasing = T),]
  
  # ---------------------------------- save cell assignments as a text file
  write.table(assignmentTable_dedup, paste0(output.dir_perSample, '/', project.id, '_', gex.library.id, '_', 'Celltype_assignment_toCells.txt'), quote = F, row.names = F, sep = '\t')
  write.table(assignmentTable_dedup_ClusterNumber, paste0(output.dir_perSample, '/', project.id, '_', gex.library.id, '_', 'Celltype_assignment_toCells_PlusClusters.txt'), quote = F, row.names = F, sep = '\t')
  
  # ---------------------------------- superimpose the cell type on cell Range clusters [tSNE]
  
  if(length(grep('grouped', input.dir))==1) # catch integer(0)
  {
    
    tSNE_Cellranger = fread(paste0(project.id,'/10X-gex-grouped/',gex.library.id,'/outs/analysis/tsne/2_components/projection.csv'), stringsAsFactors = F)
    
  }else{
    tSNE_Cellranger = fread(paste0(project.id,'/10X-gex/',gex.library.id,'/outs/analysis/tsne/2_components/projection.csv'), stringsAsFactors = F)
  }
  
  tSNE_Cellranger = as.data.frame(tSNE_Cellranger)
  row.names(tSNE_Cellranger) = gsub('-.*', '',  tSNE_Cellranger$Barcode)
  
  # ---------------------------------- keep the enriched cells [?]
  tSNE_Cellranger = tSNE_Cellranger[which(row.names(tSNE_Cellranger) %in% assignmentTable_dedup_ClusterNumber$cell),]
  dim(tSNE_Cellranger)
  
  # ----------------------------------
  cellsTsne = as.matrix(tSNE_Cellranger[,c(2,3)])
  
  # ---------------------------------- use the Cell-Ranger tSNE [test plot]
  # dev.off()
  # plot(cellsTsne, pch=16, cex=.3)
  #  
  
  # ---------------------------------- select 6 top dominant enriched cell types
  selectedThresholds <- getThresholdSelected(cells_assignment)
  selectedCellTypes = unique(assignmentTable_dedup_ClusterNumber$geneSet)[1:6]
  selectedThresholds = selectedThresholds[which(names(selectedThresholds)%in%selectedCellTypes)]
  
  # ---------------------------------- visualization of tSNE plot + functional annotation
  set.seed(123)
  pdf(paste0(output.dir_perSample, '/', project.id, '_', gex.library.id, '_', 'tSNE_Plot.pdf'), height = 8, width = 10, useDingbats = F)
  par(mfrow=c(2,3)) # Splits the plot into two rows and three columns
  
  for(geneSetName in selectedCellTypes)
  {
    nBreaks <- 5 # Number of levels in the color palettes
    # Color palette for the cells that do not pass the threshold
    colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
    # Color palette for the cells that pass the threshold
    colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
    
    passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
    
    if(sum(passThreshold) > 0 )
    {
      aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
      
      # Assign cell color
      cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                     setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
      
      plot(cellsTsne, main=geneSetName,
           sub="Pink/red cells pass the threshold",
           col=cellColor[rownames(cellsTsne)], pch=16, cex = 0.8) 
    }
  }
  dev.off()  
  
  # ---------------------------------- superimpose the cell type on cell Range clusters [UMAP]
  if(length(grep('grouped', input.dir))==1) # catch integer(0)
  {
    UMAP_Cellranger = fread(paste0(project.id,'/10X-gex-grouped/',gex.library.id,'/outs/analysis/umap/2_components/projection.csv'), stringsAsFactors = F)
    
  }else{
    UMAP_Cellranger = fread(paste0(project.id,'/10X-gex/',gex.library.id,'/outs/analysis/umap/2_components/projection.csv'), stringsAsFactors = F)
  }
  
  UMAP_Cellranger = as.data.frame(UMAP_Cellranger)
  row.names(UMAP_Cellranger) = gsub('-.*', '',  UMAP_Cellranger$Barcode)
  
  # ---------------------------------- keep the enriched cells [?]
  UMAP_Cellranger = UMAP_Cellranger[which(row.names(UMAP_Cellranger) %in% assignmentTable_dedup_ClusterNumber$cell),]
  dim(UMAP_Cellranger)
  
  # ----------------------------------
  cellsUMAP = as.matrix(UMAP_Cellranger[,c(2,3)])
  
  # ---------------------------------- select 6 top dominant enriched cell types
  selectedThresholds <- getThresholdSelected(cells_assignment)
  selectedCellTypes = unique(assignmentTable_dedup_ClusterNumber$geneSet)[1:6]
  selectedThresholds = selectedThresholds[which(names(selectedThresholds)%in%selectedCellTypes)]
  
  # ---------------------------------- visualization of UMAP plot + functional annotation
  set.seed(123)
  pdf(paste0(output.dir_perSample, '/', project.id, '_', gex.library.id, '_', 'UMAP_Plot.pdf'), height = 8, width = 10, useDingbats = F)
  par(mfrow=c(2,3)) # Splits the plot into two rows and three columns
  
  for(geneSetName in selectedCellTypes)
  {
    nBreaks <- 5 # Number of levels in the color palettes
    # Color palette for the cells that do not pass the threshold
    colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
    # Color palette for the cells that pass the threshold
    colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
    
    passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
    
    if(sum(passThreshold) > 0 )
    {
      aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
      
      # Assign cell color
      cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                     setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
      
      plot(cellsUMAP, main=geneSetName,
           sub="Pink/red cells pass the threshold",
           col=cellColor[rownames(cellsUMAP)], pch=16, cex = 0.8) 
    }
  }
  dev.off()
  
  ## -------------------------------------------------------------------------- knee plot
  matrix_dir = paste0(project_id, input.dir, '/outs/raw_feature_bc_matrix/')
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  
  # Computing barcode rank statistics:
  br.out <- barcodeRanks(mat)
  
  pdf(paste0(output.dir_perSample, '/', project.id, '_', gex.library.id, '_', 'Knee_Plot.pdf'), height = 8, width = 8, useDingbats = F)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Barcodes", ylab="UMI counts")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  dev.off()
  
  ## --------------------------------------------------------------------------
  date()
  sessionInfo()
  
  print("You can find the results at: ")
  print(paste0(getwd(),'/', project_id, input.dir))
  
}

