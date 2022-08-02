# Bulk AYA Gliloma analysis
# Last updated: 15/07/2022

########## Initialize ##########
## Install global packages ##
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tximport)

## Install custom functions
source("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/dfci-pnoc/bulkRNAseq-scripts.R")

## Set paths
wd <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/")
fwd <- c(paste0(wd, "Figures/"))

########## Load Data ##########
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

##### Load bulk RNAseq Data
cm_fc <- MakeCountMatrixFromFeatureCounts(filePath = "~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/star-fc/")
cm_rsem <- MakeCountMatrixFromRSEM(filePath = "~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/star-rsem/cm_counts_counts-rsem.csv")
txi.rsem <- MakeCountMatrixWithTximport(filePath = "~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/star-rsem/", method = "rsem")
txi.salmon <- MakeCountMatrixWithTximport(filePath = "~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/salmon/", method = "salmon")

library(DESeq2)
meta_for_dds <- meta %>% tibble::column_to_rownames(var = "SampleID.new")
sample_order <- txi.rsem$counts %>% colnames()
meta_for_dds <- meta_for_dds[sample_order,]
dds.rsem <- DESeqDataSetFromTximport(txi = txi.rsem, colData = meta_for_dds, design = ~1)
dds.salmon <- DESeqDataSetFromTximport(txi = txi.salmon, colData = meta_for_dds, design = ~1)

cm_fc <- Ensembl2Symbol(countMatrix = cm_fc) %>%
  column_to_rownames(var = "hgnc_symbol")
cm_rsem <- Ensembl2Symbol(countMatrix = cm_rsem) %>% 
  column_to_rownames(var = "hgnc_symbol")
cm_salmon <- Ensembl2Symbol(countMatrix = dds.salmon@assays@data$counts) %>% 
  column_to_rownames(var = "hgnc_symbol")

shared_genes <- intersect(rownames(cm_rsem), rownames(cm_salmon))

##### Load QC Files
filePath_qc <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/star-fc/star-fc.QC/")
qcFile_fc <- MakeQCFile(filePath.QC = filePath_qc, tool = "featureCounts")
qcFile_rsem <- MakeQCFile(filePath.QC = "~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/qc_counts-rsem.csv",
                          tool = "rsem")

##### Load gene set files
syn_genesets <- readRDS("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/General RDS/synaptic_markers_list.Rds")
source("~/Dropbox (Partners HealthCare)/Filbin lab/Shared/Scripts/Survival_HelperFunctions.R")

########## Process bulk RNA-seq ##########
##### Compare expression profiles of salmon and rsem #####
comparison <- data.frame(patient = as.character(), cor = as.numeric())
for (i in 1:length(colnames(cm_salmon))) {
  patient <- colnames(cm_salmon)[i]
  cor <- cor(cm_salmon[shared_genes, patient], cm_rsem[shared_genes, patient])
  df <- data.frame(patient, cor)
  comparison <- rbind(comparison, df)
}

salmon.rsem_Cor <- cor(cm_salmon[shared_genes,], cm_rsem[shared_genes,])
pheatmap(salmon.rsem_Cor, cluster_cols = FALSE, cluster_rows = FALSE)

densityDF <- data.frame(sum = c(colSums(cm_salmon), colSums(cm_rsem), colSums(countFile_fc)), method = c(rep("salmon", 23), rep("rsem", 23), rep("featureCounts", 23)))
ggplot(densityDF, aes(x = sum, fill = method)) +
  geom_density(alpha = 0.4) +
  xlim(c(min(densityDF$sum)-10^7/2, max(densityDF$sum)+10^7/2)) +
  xlab(label = "Total counts")
### something weird is going on with featureCounts... the alignment rate is much lower than rsem, and rsem and salmon are very comparable

cpm_rsem <- apply(cm_rsem, 2, function(x) (x/sum(x))*1000000)

### Use RSEM pipeline from here on our...

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
input <- CenterCountsMatrix(countMatrix = cpm_rsem)
meta2 <- ScoreSignatureAndSplit(countMatrix.center = input$countMatrix_center,
                                countMatrix.mean = input$countMatrix_mean,
                                metaData = meta,
                                s = "OLIG2",
                                geneSetName = "test",
                                splitBy = "Median",
                                verbose = TRUE)
test <- ScoreSurvival(metaData = meta2, geneSetName = "test")
test$plot

## Examine glutamateric gene signatures
gluta.sig <- unlist(syn_genesets$curatedSynaptic$gluta)
geneSignif <- data.frame(gene = as.character(), pvalue = as.numeric())
for (gene in gluta.sig) {
  metaGene <- ScoreSignatureAndSplit(countMatrix.center = input$countMatrix_center,
                                     countMatrix.mean = input$countMatrix_mean,
                                     metaData = meta,
                                     s = gene,
                                     geneSetName = gene,
                                     splitBy = "Median",
                                     verbose = FALSE)
  test <- ScoreSurvival(metaData = metaGene,
                        geneSetName = gene)
  pval <- as.numeric(strsplit(x = test$p_val, split = " = ")[[1]][2])
  geneSignif <- rbind(geneSignif, data.frame(gene = gene, pvalue = pval))
  
}


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



