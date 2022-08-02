

# load("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/Reference data/TCGA Jones/TCGA_Jones_idhwt_YA_clinical.Robj")
# load("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/Reference data/TCGA Jones/TCGA_Jones_idhwt_YA_survival.Robj")
# load("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/Reference data/TCGA Jones/TCGAJones_idhwt_YA_centeredCounts.Robj")

load("~/Dropbox (Partners HealthCare)/Filbin lab/Jacob/Projects/AYA Bulk RNA-seq/Reference data/TCGA_ClinicalMetadata.Robj")
meta <- clinical_list[["idh1_YA"]]

load("~/Dropbox (Partners HealthCare)/Filbin lab/Shared/PublishedDatasets/BulkRNAseq/TCGA/AllTCGACohorts_CountsFPKMTPM_TCGAIDs_SurvivalMatched.Robj")
cm <- all_cm$counts$idh1_YA_counts

shared_sample <- intersect(colnames(cm), meta$case_submitter_id)
cm <- cm[, shared_sample]
meta <- meta[meta$case_submitter_id %in% shared_sample, ] %>% distinct(case_submitter_id, .keep_all = TRUE)

length(unique(b$case_id))