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
countFile_fc <- MakeCountMatrixFromFeatureCounts(filePath = filePath_fc)

filePath_rsem <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/star-rsem/")
countFile_rsem <- MakeCountMatrixFromRSEM(filePath = filePath_rsem)

countFile_fc <- Ensembl2Symbol(count_matrix = countFile_fc)
countFile_fc <- countFile_fc %>% column_to_rownames(var = "hgnc_symbol")

countFile_rsem <- Ensembl2Symbol(count_matrix = countFile_rsem)
countFile_rsem <- countFile_rsem %>% column_to_rownames(var = "hgnc_symbol")

shared_genes <- intersect(rownames(countFile_fc), rownames(countFile_rsem))

##### Load QC Files
filePath_qc <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/star-fc/star-fc.QC/")
qcFile_fc <- MakeQCFileFromFeatureCountsOutput(filePath.QC = filePath_qc)

##### Load AYA clinical data
meta <- read_excel("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Data/AYA Bulk RNA-seq/Deidentified RNAseq Patient Clinical Data.xlsx",
                   skip = 1) 
meta <- meta %>%
  mutate(SampleID.new = paste0("s", sapply(strsplit(as.character(meta$SampleID), "-"), `[`, 1))) %>% 
  # only keeps sample ID before "-"
  filter(SampleID.new %in% colnames(countFile_fc)) %>%
  filter(SampleID != "94410-1*") %>%
  mutate(time = as.numeric(`Overall Survival (Years)`)) %>% 
  mutate(pathology16 = `Pathology at Diagnosis`) %>%
  mutate(pathology21 = ifelse(`Pathology at Diagnosis` %in% c("Anaplastic Astrocytoma, WHO Grade III, IDH-Mutant", "Diffuse Astrocytoma, WHO Grade II, IDH-Mutant"), "Astrocytoma, IDH-Mutant",
                              ifelse(`Pathology at Diagnosis` == "Oligodendroglioma, WHO Grade II", "Oligodendroglioma, IDH-Mutant", "Glioblastoma"))) %>%
  mutate(status = ifelse(`Vital Status` == "Alive", 1, 2)) # this will be used in `status` in survival package

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
  xlim(c(min-10^7/2, max+10^7/2)) +
  xlab(label = "Total gene expression")

cpm_fc <- apply(countFile_fc, 2, function(x) (x/sum(x))*1000000)
cpm_rsem <- apply(countFile_rsem, 2, function(x) (x/sum(x))*1000000)
for (sample in colnames(cpm_fc)) {
  print(sample)
  print(cor(cpm_fc[shared_genes, sample], cpm_rsem[shared_genes, sample]))
}
# Use featureCounts output from here on out
cpm <- cpm_fc

rm(cpm_rsem)
rm(countFile_rsem)
rm(filePath_rsem)

##### Qualtiy control on featureCounts data
qcPlots <- MakeQCPlotsFromQCData(qcFile = qcFile_fc, countMatrix = cm)
qcPlots[[1]]
qcPlots[[2]]
qcPlots[[3]]

########## Survival on clinical metadata ##########
library(survival)
library(survminer)

fit <- survfit(Surv(time = time, status) ~ pathology16, data = meta)
pvalue <- surv_pvalue(fit, meta)$pval.txt
ggsurvplot(fit,
           pval = TRUE,
           risk.table = TRUE, 
           risk.table.y.text = FALSE,
           risk.table.height = 0.25,
           legend.labs = c("Anaplastic Astro.", "Diffuse Astro.", "Glioblastoma", "Oligodendroglioma"))

########## Score sample by expression ##########
input <- CenterCountsMatrix(cm = cpm)
meta$geneSig <- scoreSignature(input$cm_center, input$cm_mean, s = syn_genesets$curatedSynaptic$gluta, simple = FALSE, verbose = TRUE)
ggplot(meta, aes(x = pathology16, y = geneSig)) +
  geom_boxplot() +
  geom_jitter()

meta$gaba <- scoreSignature(input$cm_center, input$cm_mean, s = syn_genesets$curatedSynaptic$gaba, simple = FALSE, verbose = TRUE)
meta$general_synaptic <- scoreSignature(input$cm_center, input$cm_mean, s = syn_genesets$curatedSynaptic$general_synaptic, simple = FALSE, verbose = TRUE)
intersect(syn_genesets$curatedSynaptic$gaba, syn_genesets$curatedSynaptic$gluta) # these gene sets do not overlap
plot(meta$gluta, meta$gaba)
hist(meta$gluta)
RunSurvivalAnalysis()



