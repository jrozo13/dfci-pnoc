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

## Set paths
wd <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/")
fwd <- c(paste0(wd, "Figures/"))

########## Load Data ##########
##### Load bulk RNAseq Data
filePath_fc <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/star-fc/")
countFile_fc <- MakeCountMatrixFromFeatureCounts(filePath = filePath_fc)

countFile_rsem <- MakeCountMatrixFromRSEM(filePath = "~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/cm_counts_counts-rsem.csv")

countFile_fc <- Ensembl2Symbol(countMatrix = countFile_fc) %>% 
  column_to_rownames(var = "hgnc_symbol")

countFile_rsem <- Ensembl2Symbol(countMatrix = countFile_rsem) %>% 
  column_to_rownames(var = "hgnc_symbol")

shared_genes <- intersect(rownames(countFile_fc), rownames(countFile_rsem))

##### Load QC Files
filePath_qc <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/star-fc/star-fc.QC/")
qcFile_fc <- MakeQCFile(filePath.QC = filePath_qc, tool = "featureCounts")
qcFile_rsem <- MakeQCFile(filePath.QC = "~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/qc_counts-rsem.csv",
                          tool = "rsem")
##### Load AYA clinical data
meta <- read_excel("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/Deidentified RNAseq Patient Clinical Data.xlsx",
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
syn_genesets <- readRDS("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/General RDS/synaptic_markers_list.Rds")
source("~/Dropbox (Partners HealthCare)/Filbin lab/Shared/Scripts/Survival_HelperFunctions.R")

########## Process bulk RNA-seq ##########
##### Compare expression profiles of featureCounts and rsem #####
comparison <- data.frame(patient = as.character(), cor = as.numeric())
for (i in 1:length(colnames(countFile_fc))) {
  patient <- colnames(countFile_fc)[i]
  cor <- cor(countFile_fc[shared_genes, patient], countFile_rsem[shared_genes, patient])
  df <- data.frame(patient, cor)
  comparison <- rbind(comparison, df)
}

fc.rsem_Cor <- cor(countFile_fc[shared_genes,], countFile_rsem[shared_genes,])
pheatmap(fc.rsem_Cor, cluster_cols = FALSE, cluster_rows = FALSE)

densityDF <- data.frame(sum = c(colSums(countFile_fc), colSums(countFile_rsem)), method = c(rep("featureCounts", 23), rep("rsem", 23)))
ggplot(densityDF, aes(x = sum, fill = method)) +
  geom_density(alpha = 0.4) +
  xlim(c(min(densityDF$sum)-10^7/2, max(densityDF$sum)+10^7/2)) +
  xlab(label = "Total counts")

cpm_fc <- apply(countFile_fc, 2, function(x) (x/sum(x))*1000000)
cpm_rsem <- apply(countFile_rsem, 2, function(x) (x/sum(x))*1000000)
comparison <- data.frame(patient = as.character(), cor = as.numeric())
for (i in 1:length(colnames(cpm_fc))) {
  patient <- colnames(cpm_fc)[i]
  cor <- cor(cpm_fc[shared_genes, patient], cpm_rsem[shared_genes, patient])
  df <- data.frame(patient, cor)
  comparison <- rbind(comparison, df)
}

##### Analysis of quality control data #####
qcPlots_rsem <- MakeQCPlots(qcFile = qcFile_rsem, countMatrix = cpm_rsem)
qcPlots_rsem.ribo <- qcPlots_rsem$pRibo
qcPlots_rsem.numAligned <- qcPlots_rsem$nAligned
qcPlots_rsem.percAligned <- qcPlots_rsem$pAligned
print(qcPlots_fc.ribo, newpage = FALSE)

pdf(paste0(fwd, "BulkRNAseqQC.rsem.pdf"), width = 15, height = 5)
ggarrange(qcPlots_rsem.ribo, qcPlots_rsem.numAligned, qcPlots_rsem.percAligned,
                                  ncol = 3, nrow = 1)
dev.off()

qcPlots_fc <- MakeQCPlots(qcFile = qcFile_fc, countMatrix = cpm_fc)
qcPlots_fc.ribo <- qcPlots_fc$pRibo
qcPlots_fc.numAligned <- qcPlots_fc$nAligned
qcPlots_fc.percAligned <- qcPlots_fc$pAligned

pdf(paste0(fwd, "BulkRNAseqQC.fc.pdf"), width = 15, height = 5)
ggarrange(qcPlots_fc.ribo, qcPlots_fc.numAligned, qcPlots_fc.percAligned,
          ncol = 3, nrow = 1)
dev.off()

########## Survival on clinical metadata ##########
library(survival)
library(survminer)

fit.who16 <- survfit(Surv(time = time, status) ~ pathology16, data = meta)
pvalue <- surv_pvalue(fit.who16, meta)$pval.txt
survivalPlot.who16 <- ggsurvplot(fit.who16,
           risk.table = TRUE, 
           risk.table.y.text = FALSE,
           risk.table.height = 0.25,
           legend.labs = c("Anaplastic Astro.", "Diffuse Astro.", "Glioblastoma", "Oligodendroglioma"))

pdf(paste0(fwd, "SurvivalCurve.WHO16Diagnosis.pdf"), width = 7, height = 7)
print(survivalPlot.who16, newpage = FALSE)
dev.off()

fit.who21 <- survfit(Surv(time = time, status) ~ pathology21, data = meta)
pvalue <- surv_pvalue(fit.who21, meta)$pval.txt
survivalPlot.who21 <- ggsurvplot(fit.who21,
                                 risk.table = TRUE, 
                                 risk.table.y.text = FALSE,
                                 risk.table.height = 0.25,
                                 legend.labs = c("Astrocytoma", "Glioblastoma", "Oligodendroglioma"))

pdf(paste0(fwd, "SurvivalCurve.WHO21Diagnosis.pdf"), width = 7, height = 7)
print(survivalPlot.who21, newpage = FALSE)
dev.off()

########## Score sample by expression ##########
meta2 <- ScoreSignatureAndSplit(countMatrix.center = input$countMatrix_center,
                                countMatrix.mean = input$countMatrix_mean,
                                metaData = meta,
                                s = syn_genesets$msigdb$glutamatergic_synapse,
                                geneSetName = "test",
                                splitBy = "Median",
                                verbose = TRUE)
test <- ScoreSurvival(metaData = meta2, geneSetName = "test")
test$plot



scoreDF <- meta %>% dplyr::select(SampleID.new, pathology16, pathology21, time, status)
input <- CenterCountsMatrix(countMatrix = cpm_rsem)
ggplot(scoreDF, aes(x = pathology16, y = geneSig)) +
  geom_boxplot() +
  geom_jitter()

scoreDF$gluta <- ScoreSignature_old(input$countMatrix_center, input$countMatrix_mean, s = syn_genesets$curatedSynaptic$gluta, simple = FALSE, verbose = TRUE)

markerList <- list(gluta = syn_genesets$curatedSynaptic$gluta)
scoreSplitDF <- SplitHighLow(scoreDF = scoreDF, signatureList = markerList, sampleIDColumn = "SampleID.new", splitBy = "Median")
a <- ScoreSurvival(cbind(scoreDF, scoreSplitDF), geneSetName = "gluta")
fit <- survfit(Surv(scoreDF$time, scoreDF$status) ~ scoreSplitDF$gluta_Cluster)
surv_pvalue(fit, data = scoreSplitDF)$pval.txt
plot1 <- ggsurvplot(fit,
                    data = metaData,
                    pval = TRUE,
                    risk.table = TRUE,
                    risk.table.y.text = FALSE,
                    risk.table.height = 0.25)
a$plot



