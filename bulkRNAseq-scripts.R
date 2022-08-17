######################################################
##### Functions useful for bulk RNA-seq analysis #####
######################################################

##### Load important packages #####
library(tidyverse)
library(readr)
library(biomaRt)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

##### Function to read count outputs from featureCounts #####
MakeCountMatrixFromFeatureCounts <- function(filePath) {
  files <- list.files(path = filePath, pattern = ".txt")
  
  for (i in 1:length(files)) {
    file <- files[i]
    print(file)
    sampleName <- sapply(strsplit(file,"_"), getElement, 1)
    sampleCounts <- read.table(file = paste0(filePath, file), header = TRUE)
    colnames(sampleCounts) <- c("Gene", sampleName)
    
    sampleCounts <- sampleCounts %>% separate("Gene", c("Gene", NA))
    
    if (i == 1){
      countFile <- sampleCounts
    } else {
      countFile <- merge(countFile, sampleCounts, by = "Gene")
    }
  }
  return(countFile)
}

##### Function to read count outputs from rsem #####
MakeCountMatrixFromRSEM <- function(filePath) {
  countFile <- read.csv(filePath)
  samples <- colnames(countFile)
  for (i in 2:length(samples)) {
    sample <- samples[i]
    sampleName <- sapply(strsplit(sample,"[.]"), getElement, 1)
    colnames(countFile)[i] <- sampleName
  }
  colnames(countFile)[1] <- "Gene"
  countFile <- countFile %>% separate("Gene", c("Gene", NA))
  return(countFile)
}

##### Function to read count outputs from salmon #####
MakeCountMatrixWithTximport <- function(filePath, method) {
  # mart <- readRDS("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/General RDS/ensembl_to_symbol.RDS")
  # tx_list <- getBM(filters = "ensembl_transcript_id", 
  #                  attributes= c("ensembl_transcript_id", "ensembl_gene_id"),
  #                  values = sampleCounts$Transcript,
  #                  mart = mart)
  # save(mart, tx_list, file = "/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/General RDS/tx_to_ensembl.RData")
  
  # load(file = "/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/General RDS/tx_to_ensembl.RData")
  load(file = "/Users/jrozowsk/Dropbox (Partners HealthCare)/Filbin/Jacob/Projects/General RDS/tx_to_ensembl.RData")
  
  if (method == "rsem") {
    files <- list.files(path = filePath, pattern = ".isoforms.results")
    filesFullPath <- paste0(filePath, files)
    
    names(filesFullPath) <- sapply(strsplit(files, "[.]"), getElement, 1)
    
    txi <- tximport(filesFullPath, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx_list, ignoreTxVersion = TRUE)
  }
  
  if (method == "salmon") {
    filePath = "~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/salmon/"
    files <- list.files(path = filePath, pattern = "quant.txt")
    filesFullPath <- paste0(filePath, files)
    
    names(filesFullPath) <- sapply(strsplit(files, "[.]"), getElement, 1)
    
    txi <- tximport(filesFullPath, type = "salmon", tx2gene = tx_list, ignoreTxVersion = TRUE)
  }
  
  countFile <- as.data.frame(txi$counts) %>% tibble::rownames_to_column(var = "Gene")
  countFile$Gene <- sapply(strsplit(countFile$Gene, "[.]"), getElement, 1)
  
  return(txi)
}

##### Function to make QC file from featureCounts.summary directory OR rsem quality file #####
## Note: if tool = featureCounts, provide directory; if tool = rsem, provide csv file
MakeQCFile <- function(filePath.QC, tool) {
  if (tool == "featureCounts") {
  files <- list.files(path = filePath.QC, pattern = ".txt.summary")
  for (i in 1:length(files)) {
    file <- files[i]
    print(file)
    sampleName <- sapply(strsplit(file,"_"), getElement, 1)
    sampleQC <- read.table(file = paste0(filePath.QC, file), header = TRUE)
    colnames(sampleQC) <- c("Metric", sampleName)
    if (i == 1){
      qcFile <- sampleQC
    } else {
      qcFile <- merge(qcFile, sampleQC, by = "Metric")
    }
  }
  qcFile <- qcFile %>% column_to_rownames(var = "Metric")
  qcFile <- data.frame(t(qcFile))

  qcFile$nAligned <- qcFile$Assigned
  qcFile$Assigned <- NULL
  
  qcFile$nTotal <- rowSums(qcFile)
  }
  
  if (tool == "rsem") {
    qcFile <- read.csv(filePath.QC)
    qcFile <- qcFile %>% separate("sample", c("sample", NA)) %>%
      column_to_rownames(var = "sample")
  }
  
  qcFile <- qcFile %>% dplyr::select(nTotal, nAligned)
  return(qcFile)
}

##### Function to make QC plots from qcFile and count matrix #####
MakeQCPlots <- function(qcFile, countMatrix) {
  riboGenes <- read_csv("~/Dropbox (Partners HealthCare)/Filbin lab/Shared/Marker_genes/RibosomalGenes.csv", col_names = "Gene", show_col_types = FALSE) %>% pull(Gene)
  nRiboCounts <- countMatrix[rownames(countMatrix) %in% riboGenes,]
  riboQC <- data.frame(RiboPercent = colSums(nRiboCounts)/colSums(countMatrix), CountSum = colSums(countMatrix))
  qcFile <- cbind(qcFile, riboQC) %>% rownames_to_column(var = "sample")
  
  plot1 <- ggplot(qcFile, aes(x = sample, y = RiboPercent))+ 
    geom_bar(position="dodge", stat="identity")+
    ylim(0, 0.25)+
    xlab("")+
    ylab("Percent of counts")+
    ggtitle("Percentage of ribosomal counts") +
    theme(axis.text=element_text(color = "black", face = "bold"),
          axis.title.y = element_text(color = "black", face = "bold", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot2 <- ggplot(qcFile, aes(x = sample, y = nAligned))+ 
    geom_bar(position="dodge", stat="identity")+
    geom_hline(yintercept = 10^7, color="red", linetype="dashed", size=1.5)+
    xlab("")+
    ylab("Number of reads")+
    ggtitle("Number of aligned reads")+
    theme(axis.text=element_text(color="black", face="bold"),
          axis.title.y = element_text(color="black", face="bold", size=12),
          axis.text.x = element_text(angle=45, hjust=1))
  
  plot3 <- ggplot(qcFile, aes(x = sample, y = nAligned/nTotal))+ 
    geom_bar(position="dodge", stat="identity")+
    xlab("")+
    ylim(c(0,1)) +
    ylab("Percentage of total reads")+
    ggtitle("Percentage of aligned reads")+
    theme(axis.text=element_text(color="black", face="bold"),
          axis.title.y = element_text(color="black", face="bold", size=12),
          axis.text.x = element_text(angle=45, hjust=1))
  
  return(list(pRibo = plot1, nAligned = plot2, pAligned = plot3))
}

##### Function to go from ENSEMBL to Gene Symbol #####
Ensembl2Symbol <- function(countMatrix) {
  # mart was created: July 14, 2022
  # mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  mart <- readRDS("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/General RDS/ensembl_to_symbol.RDS")
  # countMatrix <- cm_rsem
  colnames(countMatrix)[1] <- "ensembl_gene_id"
  G_list <- getBM(filters = "ensembl_gene_id", 
                  attributes= c("ensembl_gene_id", "hgnc_symbol"),
                  values = countMatrix$ensembl_gene_id,
                  mart = mart)
  countMatrix <- merge(G_list, countMatrix, by = "ensembl_gene_id")
  countMatrix <- countMatrix %>% 
    subset(select = -c(ensembl_gene_id)) %>%
    filter(hgnc_symbol != "")
  
  countMatrix <- aggregate(. ~ hgnc_symbol, data = countMatrix, FUN = sum)
  
  return(countMatrix)
}

##### Function to make input count matrix for ScoreSignature #####
## Normalize and center count matrix, for input into ScoreSignature
NormCenter<-function(cm, scale_factor=10, log_base=2){
  cm = as.matrix(cm)
  cm_norm = log(cm/scale_factor+1, base=log_base)
  cm_norm_center = cm_norm-rowMeans(cm_norm)
  result = list()
  result[["raw_data"]] = cm
  result[["norm_data"]] = cm_norm
  result[["center_data"]] = cm_norm_center
  return (result)
}

CenterCountsMatrix <- function(countMatrix, alreadyCentered = FALSE){
  countMatrix_list <- list()
  if(alreadyCentered) {
    countMatrix_list$tpm <- NA
    countMatrix_list$countMatrix_center<- countMatrix
    countMatrix_list$countMatrix_mean<- rowMeans(countMatrix_list$countMatrix_center)
  } else {
    countMatrix_list$tpm <- countMatrix
    normCenter <- NormCenter(countMatrix)
    countMatrix_list$countMatrix_center<- normCenter$center_data
    countMatrix_list$countMatrix_mean<- rowMeans(log2(countMatrix_list$tpm + 1))
  }
  return(countMatrix_list)
}

##### Function to score bulk RNA-seq samples with gene-set and split for survival analysis #####
## countMatrix.center: centered relative expression (output of CenterCountsMatrix)
## countMatrix.mean: average of relative expression of each gene (log2 transformed; output of CenterCountsMatrix)
## s: gene signature in vector form
## geneSetName: what to name gene signature
## metaData: metaData file where rownames are sample names/IDs
## n: number of genes with closest average expression for control genesets, default = 100 
## splitBy: can split by Median (default) or Quartiles (type "Quartile")
## simple: whether use average, default  = FALSE
ScoreSignatureAndSplit <- function(countMatrix.center, countMatrix.mean, s, geneSetName, metaData,
                                   n = 100, splitBy = "Median", simple = FALSE, verbose = FALSE) {
  if (!splitBy %in% c("Median", "Quantile")) {
    stop("splitBy parameter not in: 'Median' or 'Quantile'")
  }
  
  if(verbose) {
    message("cells: ", ncol(countMatrix.center))
    message("genes: ", nrow(countMatrix.center))
    message("genes in signature: ", length(s))
    message("Using simple average?", simple)
    message("processing...")
  }
  
  s <- intersect(rownames(countMatrix.center), s)
  message("genes in signature, and also in this dataset: ", length(s))
  
  if (simple){
    sigScore <- colMeans(v[s,])
  } else {
    sigScore <- colMeans(do.call(rbind, lapply(s, function(g) {
      if(verbose) message(".", appendLF = FALSE)
      g.n <- names(sort(abs(countMatrix.mean[g] - countMatrix.mean))[2:(n+1)])
      countMatrix.center[g, ] - colMeans(countMatrix.center[g.n, ])
    })))
  }
  
  if(verbose) message(" done")
  
  sigScore <- sigScore %>% as.data.frame()
  colnames(sigScore)[1] <- "geneSet"
  
  # if (splitBy == "Quartile") { n = 4 }
  # if (splitBy == "Median") { n = 2 }
  # 
  # scoreSplit <- sigScore %>% mutate(clusters = ntile(sigScore[[geneSetName]], n))
  # 
  if(splitBy == "Quantile"){
    params <- quantile(sigScore$geneSet)
    scoreSplit <- sigScore %>% mutate(clusters = ifelse(geneSet > params[["75%"]], "High",
                                                        ifelse(geneSet < params[["25%"]], "Low", "Mid")))
    }
  ## Split based on median
  if(splitBy == "Median"){
    params <- median(sigScore$geneSet)
    scoreSplit <- sigScore %>% mutate(clusters = ifelse(geneSet < params, "Low", "High"))
  }
  
  colnames(scoreSplit) <- c(geneSetName, paste0(geneSetName, "_Cluster"))
  scoreSplit <- scoreSplit %>% rownames_to_column(var = "SampleID.new")
  
  metaData <- merge(metaData, scoreSplit, by = "SampleID.new")
  
  return(metaData)
}

ScoreSurvival <- function(metaData, geneSetName){
  ## Prep data for survival analysis- need time to death, patient, vital status, and group (high/low)
 # metaData$status <- as.logical(metaData$status)
  metaData$time <- as.numeric(metaData$time)
  metaData$Marker <- metaData[,paste0(geneSetName, "_Cluster")]
  
  ## Run survival analysis
  fit <- survfit(Surv(time, status == 2) ~Marker, data = metaData)
  pvalue <- surv_pvalue(fit, metaData)$pval.txt
  paste(fit)
  plot1 <- ggsurvplot(fit,
                      data = metaData,
                      pval = TRUE,
                      risk.table = TRUE,
                      risk.table.y.text = FALSE,
                      risk.table.height = 0.25)
  return(list(plot = plot1, p_val = pvalue, metaData = metaData))
}

#########################
##### Old Functions #####
#########################

##### Function to score bulk RNA-seq samples with gene-set #####
## countMatrix.center: centered relative expression (output of CenterCountsMatrix)
## countMatrix.mean: average of relative expression of each gene (log2 transformed; output of CenterCountsMatrix)
## n: number of genes with closest average expression for control genesets, default = 100 
## simple: whether use average, default  = FALSE
# ScoreSignature_old <- function(countMatrix.center, countMatrix.mean, s, n = 100, simple = FALSE, verbose=FALSE) {
#   if(verbose) {
#     message("cells: ", ncol(countMatrix.center))
#     message("genes: ", nrow(countMatrix.center))
#     message("genes in signature: ", length(s))
#     message("Using simple average?", simple)
#     message("processing...")
#   }
# 
#   s <- intersect(rownames(countMatrix.center), s)
#   message("genes in signature, and also in this dataset: ", length(s))
# 
#   if (simple){
#     s.score <- colMeans(v[s,])
#   } else {
#     s.score <- colMeans(do.call(rbind, lapply(s, function(g) {
#       if(verbose) message(".", appendLF = FALSE)
#       g.n <- names(sort(abs(countMatrix.mean[g] - countMatrix.mean))[2:(n+1)])
#       countMatrix.center[g, ] - colMeans(countMatrix.center[g.n, ])
#     })))
#   }
# 
#   if(verbose) message(" done")
#   return(s.score)
# }

## From a dataframe of scores for each geneset, split samples into high/low expressors
## For each module score, split samples into "high" or "low" expression
## scores: a dataframe with genesets as columns, samples as rows 
## splitBy: denotes how to split samples into high/low geneset scorers
# SplitHighLow <- function(scoreDF, sampleIDColumn, signatureList, splitBy = "Median"){
#   AllClusters <- data.frame(scoreDF[[sampleIDColumn]])
#   colnames(AllClusters) <- sampleIDColumn
#   for (m in names(signatureList)){
#     print(m) # m is a geneset
#     df <- scoreDF %>% dplyr::select(sampleIDColumn, m) %>% as.data.frame()
#     ClusterName <- paste0(m, "_Cluster")
#     df[[m]] <- as.numeric(df[[m]])
# 
#     ## Alternative methods for splitting into high/low
#     ## Split into top25%/mid50%/bottom25%
#     if(splitBy == "SplitInto3Quartiles"){
#       df[,"ClusterName"] <- df[[m]] > quantile(df[[m]])["75%"]
#       high <- df[df[[m]] > quantile(df[[m]])["75%"],]; high[,"ClusterName"]<-"High"
#       low <- df[df[[m]] < quantile(df[[m]])["25%"],]; low[,"ClusterName"]<-"Low"
#       mid <- df[!(df[[m]] %in% c(high[[m]], low[[m]])),]; mid[,"ClusterName"]<-"Mid"
#       df <- rbind(high, low); df <- rbind(df, mid)
#     }
#     ## Split based on median
#     if(splitBy=="Median"){
#       df[,"ClusterName"]<- df[[m]] > median(df[[m]])
#       df[,"ClusterName"]<- gsub("FALSE", "Low",gsub("TRUE", "High", df[,"ClusterName"]))
#     }
# 
#     ## Convert to factor, with "Low" as comparison
#     #df$ClusterName<- factor(df$ClusterName, levels=c("High", "Low"))
#     colnames(df)<- c(sampleIDColumn, m, ClusterName)
#     AllClusters <- merge(AllClusters, df, by = sampleIDColumn)
#   }
#   return(AllClusters)
# }
