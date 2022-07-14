# dfci-pnoc
scripts for filbin/pnoc research projects

ensembl_to_symbol.R: this script contains functions taht convert between ensembl id and gene symbol using biomart R package. only input is a counts_matrix where the first column is ensembl id. this function adds gene symbol column and function of the gene/protein.