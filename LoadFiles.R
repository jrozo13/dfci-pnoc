# Bulk AYA Gliloma analysis: make datasets
# Last updated: 13/07/2022

########## Initialize ##########
## Install global packages ##
library(tidyverse)
library(org.Hs.eg.db)

wd <- 
setwd(wd)

########## Load Rozowsky data ##########
aya_bulkFiles <- list.files(path = "/Users/filbinlab/Documents/Jacob/Data/AYA Bulk RNA-seq/")
filePath <- c("/Users/filbinlab/Documents/Jacob/Data/AYA Bulk RNA-seq/")


for (i in 1:length(aya_bulkFiles)) {
  file <- aya_bulkFiles[i]
  print(file)
  sampleName <- sapply(strsplit(file,"_"), getElement, 1)
  sampleCounts <- read.table(file = paste0(filePath, file), header = TRUE)
  colnames(sampleCounts) <- c("Gene", sampleName)
  
  if (i == 1){
    countFile <- sampleCounts
  } else {
    countFile <- merge(countFile, sampleCounts, by = "Gene")
  }
}

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- countFile$Gene %>% data.frame()
colnames(genes) <- "Gene"
genes <- genes %>% separate("Gene", c("Gene", NA)) %>% pull(Gene)

G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id", "hgnc_symbol", "description"),
                values = genes,
                mart = mart)

countFile <- merge(countFile, G_list %>% select(ensembl_gene_id, ))


library("EnsDb.Hsapiens.v86")
symbols <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", column="SYMBOL")
