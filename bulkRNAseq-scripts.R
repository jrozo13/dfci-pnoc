##### Function to read count outputs from featureCounts #####
MakeCountMatrixFromFeatureCounts <- function(filePath) {
  library(tidyverse)
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
  library(tidyverse)
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

##### Function to make QC file from featureCounts.summary outputs #####
MakeQCFileFromFeatureCountsOutput <- function(filePath.QC) {
  library(tidyverse)
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
  qcFile$Unassigned_All <- rowSums(qcFile) - qcFile$Assigned
  return(qcFile)
}

##### Function to make QC plots from qcFile and count matrix #####
MakeQCPlotsFromQCData <- function(qcFile, countMatrix) {
  library(readr)
  riboGenes <- read_csv("~/Dropbox (Partners HealthCare)/Filbin lab/Shared/Marker_genes/RibosomalGenes.csv", col_names = "Gene", show_col_types = FALSE) %>% pull(Gene)
  nRiboCounts <- countMatrix[rownames(countMatrix) %in% riboGenes,]
  riboQC <- data.frame(RiboPercent = colSums(nRiboCounts)/colSums(countMatrix), CountSum = colSums(countMatrix))
  qcFile <- cbind(qcFile, riboQC) %>% rownames_to_column(var = "sample")
  
  plot1 <- ggplot(qcFile, aes(x=sample, y=RiboPercent))+ 
    geom_bar(position="dodge", stat="identity")+
    ylim(0, 0.25)+
    xlab("")+
    ylab("Percent of counts")+
    ggtitle("Percentage of ribosomal counts") +
    theme(axis.text=element_text(color = "black", face = "bold"),
          axis.title.y = element_text(color = "black", face = "bold", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot2 <- ggplot(qcFile, aes(x = sample, y = Assigned))+ 
    geom_bar(position="dodge", stat="identity")+
    geom_hline(yintercept = 10^7, color="red", linetype="dashed", size=1.5)+
    xlab("")+
    ylab("Number of assigned reads")+
    ggtitle("Number of reads")+
    theme(axis.text=element_text(color="black", face="bold"),
          axis.title.y = element_text(color="black", face="bold", size=12),
          axis.text.x = element_text(angle=45, hjust=1))
  
  plot3 <- ggplot(qcFile, aes(x = sample, y = Unassigned_All/(Unassigned_All + Assigned)))+ 
    geom_bar(position="dodge", stat="identity")+
    xlab("")+
    ylim(c(0,1)) +
    ylab("Percentage of total reads")+
    ggtitle("Percentage of unassigned reads")+
    theme(axis.text=element_text(color="black", face="bold"),
          axis.title.y = element_text(color="black", face="bold", size=12),
          axis.text.x = element_text(angle=45, hjust=1))
  
  return(list(plot1, plot2, plot3))
}

##### Function to go from ENSEMBL to Gene Symbol #####
Ensembl2Symbol <- function(count_matrix) {
  library(biomaRt)
  library(tidyverse)
  
  # mart was created: July 14, 2022
  # mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  mart <- readRDS("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/General RDS/ensembl_to_symbol.RDS")
  
  colnames(count_matrix)[1] <- "ensembl_gene_id"
  G_list <- getBM(filters = "ensembl_gene_id", 
                  attributes= c("ensembl_gene_id", "hgnc_symbol"),
                  values = count_matrix$ensembl_gene_id,
                  mart = mart)
  count_matrix <- merge(G_list, count_matrix, by = "ensembl_gene_id")
  count_matrix <- count_matrix %>% 
    subset(select = -c(ensembl_gene_id)) %>%
    filter(hgnc_symbol != "")
  count_matrix <- aggregate(. ~ hgnc_symbol, data = count_matrix, FUN = sum)
  
  return(count_matrix)
}

##### Function to make input count matrix for ScoreSignature #####
CenterCountsMatrix<- function(cm, alreadyCentered=FALSE){
  cm_list<- list()
  if(alreadyCentered){
    cm_list$tpm<- NA
    cm_list$cm_center<- cm
    cm_list$cm_mean<- rowMeans(cm_list$cm_center)
  } else{
    normCenter<- NormCenter(cm)
    cm_list$tpm<- cm
    cm_list$cm_center<- normCenter$center_data
    cm_list$cm_mean<- rowMeans(log2(cm_list$tpm + 1))
  }
  return(cm_list)
}

##### Function to score bulk RNA-seq samples with gene-set #####
## @param X.center centered relative expression
## @param X.mean average of relative expression of each gene (log2 transformed)
## @param n number of genes with closest average expression for control genesets, default = 100 
## @param simple whether use average, default  = FALSE
scoreSignature <- function(cm.center, cm.mean, s, n=100, simple = FALSE, verbose=FALSE) {
  if(verbose) {
    message("cells: ", ncol(cm.center))
    message("genes: ", nrow(cm.center))
    message("genes in signature: ", length(s))
    message("Using simple average?", simple)
    message("processing...")
  }
  
  s <- intersect(rownames(cm.center), s)
  message("genes in signature, and also in this dataset: ", length(s))

  if (simple){
    s.score <- colMeans(v[s,])
  } else {
    s.score <- colMeans(do.call(rbind, lapply(s, function(g) {
      if(verbose) message(".", appendLF = FALSE)
      g.n <- names(sort(abs(cm.mean[g] - cm.mean))[2:(n+1)])
      cm.center[g, ] - colMeans(cm.center[g.n, ])
    })))
  }
  
  if(verbose) message(" done")
  return(s.score)
}

## From a dataframe of scores for each geneset, split samples into high/low expressors
## For each module score, split samples into "high" or "low" expression
## Input: scores = a dataframe with genesets as columns, samples as rows 
##        splitBy = denotes how to split samples into high/low geneset scorers
SplitHighLow<- function(scores, splitBy="Median"){
  AllClusters<-data.frame(submitter_id=Scores$submitter_id)
  for (m in names(marker_list_InCohort)){
    print(m)
    df<- cbind(Scores[,"submitter_id"], Scores[,m]); df<-as.data.frame(df)
    ClusterName<- paste0(m, "_Cluster")
    df$V2<-as.numeric(df$V2)
    
    ## Alternative methods for splitting into high/low
    ## Split into top25%/mid50%/bottom25%
    if(SplitBasedOn=="SplitInto3Quartiles"){
      df[,"ClusterName"]<- df$V2>quantile(df$V2)["75%"]
      high<- df[df$V2>quantile(df$V2)["75%"],]; high[,"ClusterName"]<-"High"
      low<- df[df$V2<quantile(df$V2)["25%"],]; low[,"ClusterName"]<-"Low"
      mid<-df[!(df$V1 %in% c(high$V1, low$V1)),]; mid[,"ClusterName"]<-"Mid"
      df<-rbind(high,low); df<-rbind(df, mid)
    } 
    ## Split based on median
    if(SplitBasedOn=="Median"){
      df[,"ClusterName"]<- df$V2>median(df$V2)
      df[,"ClusterName"]<- gsub("FALSE", "Low", 
                                gsub("TRUE", "High", df[,"ClusterName"]))
    } 
    
    ## Convert to factor, with "Low" as comparison
    #df$ClusterName<- factor(df$ClusterName, levels=c("High", "Low"))
    colnames(df)<- c("submitter_id", m, ClusterName)
    AllClusters<-merge(AllClusters,df, by="submitter_id")
  }
}
