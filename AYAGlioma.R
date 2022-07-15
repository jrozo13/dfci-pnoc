# Bulk AYA Gliloma analysis
# Last updated: 15/07/2022

########## Initialize ##########
## Install global packages ##
library(tidyverse)
source("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/dfci-pnoc/bulkRNAseq-scripts.R")

########## Load AYA data ##########
filePath_fc <- c("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/star-fc/")
countFile_fc <- featurec_to_countmatrix(filePath = filePath_fc)

filePath_rsem <- c("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/star-rsem/")
countFile_rsem <- rsem_to_countmatrix(filePath = filePath_rsem)

### Convert from ENSEMBL to Gene Symbol ###
countFile_fc <- ensembl_to_symbol(count_matrix = countFile_fc)
countFile_rsem <- ensembl_to_symbol(count_matrix = countFile_rsem)