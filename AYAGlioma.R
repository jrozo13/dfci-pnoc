# Bulk AYA Gliloma analysis
# Last updated: 15/07/2022

########## Initialize ##########
## Install global packages ##
library(tidyverse)
source("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/dfci-pnoc/bulkRNAseq-scripts.R")

########## Load AYA seq data ##########
filePath_fc <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/star-fc/")
countFile_fc <- featurec_to_countmatrix(filePath = filePath_fc)

filePath_rsem <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/star-rsem/")
countFile_rsem <- rsem_to_countmatrix(filePath = filePath_rsem)

### Convert from ENSEMBL to Gene Symbol ###
countFile_fc <- ensembl_to_symbol(count_matrix = countFile_fc)
countFile_fc <- countFile_fc %>% column_to_rownames(var = "hgnc_symbol")

countFile_rsem <- ensembl_to_symbol(count_matrix = countFile_rsem)
countFile_rsem <- countFile_rsem %>% column_to_rownames(var = "hgnc_symbol")

shared_genes <- intersect(rownames(countFile_fc), rownames(countFile_rsem))

### Compare expression profiles of featureCounts and rsem ###
for (sample in colnames(countFile_fc)) {
  print(sample)
  print(cor(countFile_fc[shared_genes, sample], countFile_rsem[shared_genes, sample]))
}

library(pheatmap)
fc.rsem_Cor <- cor(countFile_fc[shared_genes,], countFile_rsem[shared_genes,])
pheatmap(fc.rsem_Cor, cluster_rows = F, cluster_cols = F)

library(ggplot2)
densityDF <- data.frame(sum = c(colSums(countFile_fc), colSums(countFile_rsem)), method = c(rep("featureCounts", 23), rep("rsem", 23)))
ggplot(densityDF, aes(x = sum, fill = method)) +
  geom_density(alpha = 0.4) +
  xlim(min-10^7/2, max+10^7/2)

########## Load AYA cliical data ##########
library(readxl)
aya_clinical <- read_excel("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/Deidentified RNAseq Patient Clinical Data.xlsx",
                           skip = 1) 
aya_clinical <- aya_clinical %>%
  mutate(SampleID.new = paste0("s", sapply(strsplit(as.character(aya_clinical$SampleID), "-"), `[`, 1))) %>% 
  # only keeps sample ID before "-"
  filter(SampleID.new %in% colnames(countFile_fc)) %>%
  filter(SampleID != "94410-1*")
# check
aya_clinical$SampleID.new %in% colnames(countFile_fc) %>% table()
colnames(countFile_fc) %in% aya_clinical$SampleID.new %>% table()


