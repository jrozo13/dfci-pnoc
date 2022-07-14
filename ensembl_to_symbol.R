##### Function to go from ENSEMBL to Gene Symbol #####
ensembl_to_symbol <- function(count_matrix) {
  library(biomaRt)
  library(tidyverse)
  
  # mart was created: July 14, 2022
  # mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  mart <- readRDS("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/General RDS/ensembl_to_symbol.RDS")
  
  colnames(count_matrix)[1] <- "ensembl_gene_id"
  G_list <- getBM(filters= "ensembl_gene_id", 
                  attributes= c("ensembl_gene_id", "hgnc_symbol", "description"),
                  values = count_matrix$ensembl_gene_id,
                  mart = mart)
  count_matrix <- merge(G_list, count_matrix, by = "ensembl_gene_id")
  return(count_matrix)
}

