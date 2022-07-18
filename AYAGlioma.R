# Bulk AYA Gliloma analysis
# Last updated: 15/07/2022

########## Initialize ##########
## Install global packages ##
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)

## Install custom functions
source("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/dfci-pnoc/bulkRNAseq-scripts.R")


########## Load Data ##########
##### Load bulk RNAseq Data
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

##### Load AYA clinical data
meta <- read_excel("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/Deidentified RNAseq Patient Clinical Data.xlsx",
                           skip = 1) 
meta <- meta %>%
  mutate(SampleID.new = paste0("s", sapply(strsplit(as.character(meta$SampleID), "-"), `[`, 1))) %>% 
  # only keeps sample ID before "-"
  filter(SampleID.new %in% colnames(countFile_fc)) %>%
  filter(SampleID != "94410-1*")
# check
meta$SampleID.new %in% colnames(countFile_fc) %>% table()
colnames(countFile_fc) %in% meta$SampleID.new %>% table()

##### Load gene set files
syn_genesets <- readRDS("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/synaptic_markers_list.Rds")
source("~/Dropbox (Partners HealthCare)/Filbin lab/Shared/Scripts/Survival_HelperFunctions.R")

########## Process bulk RNA-seq ##########
### Compare expression profiles of featureCounts and rsem ###
for (sample in colnames(countFile_fc)) {
  print(sample)
  print(cor(countFile_fc[shared_genes, sample], countFile_rsem[shared_genes, sample]))
}

fc.rsem_Cor <- cor(countFile_fc[shared_genes,], countFile_rsem[shared_genes,])
pheatmap(fc.rsem_Cor, cluster_rows = F, cluster_cols = F)

densityDF <- data.frame(sum = c(colSums(countFile_fc), colSums(countFile_rsem)), method = c(rep("featureCounts", 23), rep("rsem", 23)))
ggplot(densityDF, aes(x = sum, fill = method)) +
  geom_density(alpha = 0.4) +
  xlim(min-10^7/2, max+10^7/2) +
  xlab(label = "Total gene expression")

cpm_fc <- apply(countFile_fc, 2, function(x) (x/sum(x))*1000000)
cpm_rsem <- apply(countFile_rsem, 2, function(x) (x/sum(x))*1000000)
for (sample in colnames(cpm_fc)) {
  print(sample)
  print(cor(cpm_fc[shared_genes, sample], cpm_rsem[shared_genes, sample]))
}

########## Score sample by expression ##########
input <- NormCenter(cm = cpm_fc)
input.mean <- rowMeans(log2(input$raw_data+1))
meta$gluta <- scoreSignature(input$center_data, input.mean, s = syn_genesets$curatedSynaptic$gluta, simple = FALSE, verbose = TRUE)

RunSurvivalAnalysis()



