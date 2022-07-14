# dfci-pnoc
scripts for filbin/pnoc research projects

ensembl_to_symbol.R: this script will convert from ensembl id to gene symbol using biomart R package. only input is a counts_matrix where the first column is ensembl id. this funciton adds gene symbol column and function of the gene/protein.