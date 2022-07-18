# dfci-pnoc
scripts for filbin/pnoc research projects

bulkRNAseq-scripts.R: this file contains functions that are useful for bulkRNAseq analysis.

ensembl_to_symbol: convert between ensembl id and gene symbol using biomart R package. only input is a counts_matrix where the first column is ensembl id. this function adds gene symbol column and function of the gene/protein.

featurec_to_countmatrix: takes output files from featureCounts and makes a count matrix. only input is the directory where the .txt files are contained.

rsem_to_countmatrix: same as featurec_to_countmatrix, but for rsem output.