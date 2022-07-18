##### Function to read count outputs from featureCounts #####
featurec_to_countmatrix <- function(filePath) {
  library(tidyverse)
  files <- list.files(path = filePath)
  
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
rsem_to_countmatrix <- function(filePath) {
  library(tidyverse)
  files <- list.files(path = filePath)
  
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


##### Function to go from ENSEMBL to Gene Symbol #####
ensembl_to_symbol <- function(count_matrix) {
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

##### Function to score bulk RNA-seq samples with gene-set #####
## @param X.center centered relative expression
## @param X.mean average of relative expression of each gene (log2 transformed)
## @param n number of genes with closest average expression for control genesets, default = 100 
## @param simple whether use average, default  = FALSE
scoreSignature <- function(X.center, X.mean, s, n=100, simple = FALSE, verbose=FALSE) {
  if(verbose) {
    message("cells: ", ncol(X.center))
    message("genes: ", nrow(X.center))
    message("genes in signature: ", length(s))
    message("Using simple average?", simple)
    message("processing...")
  }
  
  s <- intersect(rownames(X.center), s)
  message("genes in signature, and also in this dataset: ", length(s))
  ##message("These genes are: ", s)
  
  if (simple){
    s.score <- colMeans(X.center[s,])
  }else{
    s.score <- colMeans(do.call(rbind, lapply(s, function(g) {
      # g <- s[2]
      # message(g)
      if(verbose) message(".", appendLF = FALSE)
      g.n <- names(sort(abs(X.mean[g] - X.mean))[2:(n+1)])
      X.center[g, ] - colMeans(X.center[g.n, ])
    })))
  }
  
  if(verbose) message(" done")
  return(s.score)
}


