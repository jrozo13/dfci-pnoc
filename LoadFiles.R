# Bulk AYA Gliloma analysis: make datasets
# Last updated: 13/07/2022

########## Initialize ##########
## Install global packages ##
library(tidyverse)

########## Load Rozowsky data ##########
filePath <- c("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/")
aya_bulkFiles <- list.files(path = filePath)

for (i in 1:length(aya_bulkFiles)) {
  file <- aya_bulkFiles[i]
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

### Convert from ENSEMBL to Gene Symbol ###
source("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/dfci-pnoc/ensembl_to_symbol.R")
countFile <- ensembl_to_symbol(count_matrix = countFile)