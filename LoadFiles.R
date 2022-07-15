# Bulk AYA Gliloma analysis: make datasets
# Last updated: 13/07/2022

########## Initialize ##########
## Install global packages ##
library(tidyverse)

########## Load Rozowsky data ##########
filePath <- c("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/star-fc/")
countFile_fc <- featurec_to_countmatrix(filePath = filePath)

### Convert from ENSEMBL to Gene Symbol ###
source("/Users/filbinlab/Dropbox (Partners HealthCare)/Filbin lab/Jacob/dfci-pnoc/ensembl_to_symbol.R")
countFile <- ensembl_to_symbol(count_matrix = countFile)