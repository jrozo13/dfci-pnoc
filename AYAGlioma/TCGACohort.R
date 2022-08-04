# Bulk AYA Gliloma analysis on TCGA Data
# Last updated: 02/08/2022

########## Initialize ##########
## Install global packages ##
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)

## Install custom functions
source("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/dfci-pnoc/bulkRNAseq-scripts.R")

## Set paths
wd <- c("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/")
fwd <- c(paste0(wd, "Figures/"))

########## Load Data ##########
##### Load AYA clinical data
load("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/Reference data/TCGA_ClinicalMetadata.Robj")
meta <- clinical_list[["idh1_YA"]]
keep_cols <- c("case_id", "case_submitter_id", "project_id", "age_at_index", "gender", "race",
               "vital_status", "primary_diagnosis", "days_to_death", "days_to_last_follow_up",
               "site_of_resection_or_biopsy", "treatment_or_therapy", "treatment_type")

##### Load bulk RNAseq Data
load("~/Dropbox (Partners HealthCare)/Filbin lab/Shared/PublishedDatasets/BulkRNAseq/TCGA/AllTCGACohorts_CountsFPKMTPM_TCGAIDs_SurvivalMatched.Robj")
cm <- all_cm$counts$idh1_YA_counts; rm(all_cm)
shared_sample <- intersect(colnames(cm), meta$case_submitter_id)

## Update metadata and count matrix with shared samples
cm <- cm[, shared_sample]
meta <- meta[meta$case_submitter_id %in% shared_sample, ] %>% 
  distinct(case_submitter_id, .keep_all = TRUE) %>% 
  dplyr::select(keep_cols) %>%
  mutate(status = ifelse(vital_status == "Alive", 1, ifelse(vital_status == "Dead", 2, 1))) %>% # this will be used in `status` in survival package
  mutate(time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death)) %>%
  mutate(SampleID.new = case_submitter_id)

##### Load gene set files
syn_genesets <- readRDS("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/General RDS/synaptic_markers_list.Rds")
# source("~/Dropbox (Partners HealthCare)/Filbin lab/Shared/Scripts/Survival_HelperFunctions.R")

########## Survival on clinical metadata ##########
library(survival)
library(survminer)

fit <- survfit(Surv(time = time, status) ~ primary_diagnosis, data = meta)
pvalue <- surv_pvalue(fit, meta)$pval.txt
survivalPlot <- ggsurvplot(fit,
                           risk.table = TRUE, 
                           risk.table.y.text = FALSE,
                           risk.table.height = 0.25)
print(survivalPlot)

pdf(paste0(fwd, "SurvivalCurve.TCGA.Diagnosis.pdf"), width = 7, height = 7)
print(survivalPlot, newpage = FALSE)
dev.off()

########## Score sample by expression ##########
input <- CenterCountsMatrix(countMatrix = cm)
## Examine glutamateric gene signatures
gluta.sig <- unlist(syn_genesets$curatedSynaptic$gluta)
geneSignif <- data.frame(gene = as.character(), pvalue = as.numeric())
meta.temp <- meta
for (gene in gluta.sig) {
  metaGene <- ScoreSignatureAndSplit(countMatrix.center = input$countMatrix_center,
                                     countMatrix.mean = input$countMatrix_mean,
                                     metaData = meta.temp,
                                     s = gene,
                                     geneSetName = gene,
                                     splitBy = "Quantile",
                                     verbose = FALSE)
  test <- ScoreSurvival(metaData = metaGene,
                        geneSetName = gene)
  pval <- as.numeric(strsplit(x = test$p_val, split = " = ")[[1]][2])
  geneSignif <- rbind(geneSignif, data.frame(gene = gene, pvalue = pval))
}

# GRIK3 and GRIK5 high expression: confers worse survival
# So goes HOX genes... why...
# These are two genes related to glutamatergic signaling

meta <- ScoreSignatureAndSplit(countMatrix.center = input$countMatrix_center,
                                countMatrix.mean = input$countMatrix_mean,
                                metaData = meta,
                                s = syn_genesets$curatedSynaptic$gluta,
                                geneSetName = "gluta",
                                splitBy = "Median",
                                verbose = TRUE)
meta <- ScoreSignatureAndSplit(countMatrix.center = input$countMatrix_center,
                               countMatrix.mean = input$countMatrix_mean,
                               metaData = meta,
                               s = syn_genesets$curatedSynaptic$gaba,
                               geneSetName = "gaba",
                               splitBy = "Median",
                               verbose = TRUE)
test <- ScoreSurvival(metaData = meta, geneSetName = "gluta")
test$plot

########## Create gene signature on GRIK3 and GRIK5 expression ##########
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

dds <- DESeqDataSetFromMatrix(countData = cm, colData = meta, design = ~IDH1YA.gs_Cluster)
dds <- DESeq(dds)
res <- na.omit(results(dds, contrast=c("IDH1YA.gs_Cluster", "High", "Low")))
res <- as.data.frame(res)
res <- res[order(res$padj),]
res$Group <- ifelse(res$log2FoldChange > 0, "High", "Low")

EnhancedVolcano(res,
                lab=rownames(res),
                x="log2FoldChange",
                y="padj",
                title="High v Low Signature Expression",
                subtitle = "",
                pCutoff=1e-10,
                FCcutoff = .75,
                col=c('black', 'black', 'black', 'red3'),
                ylab = "-Log10(p.adj)")

high.genes <- filter(res, padj < 0.05 & log2FoldChange > 0) %>% arrange(-log2FoldChange) %>% rownames() %>% unique()
low.genes <- filter(res, padj < 0.05 & log2FoldChange < 0) %>% arrange(log2FoldChange) %>% rownames() %>% unique()

gene_list <- high.genes
enrichGO(gene = gene_list,
         universe = as.character(rownames(res)),
         OrgDb = org.Hs.eg.db,
         keyType = 'SYMBOL',
         ont = "BP",
         pAdjustMethod = "BH",
         qvalueCutoff = 0.05,
         readable = FALSE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(showCategory=30)

test <- highGrik_GO@result %>%
  filter(p.adjust < 0.05)


UE_GO$Location_Enrich <- paste0(location, ".down")
assign(paste0(location, ".down"), UE_GO)



norm <- as.data.frame(apply(cm, 2, function(x){x/sum(x)*1000000}))

n_genes <- 20
top_genes <- res %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "Gene") %>% 
  group_by(Group) %>% top_n(n=n_genes, wt=padj)
top_genes <- top_genes[order(top_genes$Group),] %>% pull(Gene)



