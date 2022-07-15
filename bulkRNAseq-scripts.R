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
                  attributes= c("ensembl_gene_id", "hgnc_symbol", "description"),
                  values = count_matrix$ensembl_gene_id,
                  mart = mart)
  count_matrix <- merge(G_list, count_matrix, by = "ensembl_gene_id")
  count_matrix <- count_matrix %>% 
    subset(select = -c(ensembl_gene_id, description)) %>%
    filter(hgnc_symbol != "")
  
  return(count_matrix)
}




