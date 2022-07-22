ScoreSignature <- function(countMatrix.center, countMatrix.mean, s, geneSetName, metaData,
                           n = 100, splitBy = "Median", simple = FALSE, verbose = FALSE) {
  if(verbose) {
    message("cells: ", ncol(countMatrix.center))
    message("genes: ", nrow(countMatrix.center))
    message("genes in signature: ", length(s))
    message("Using simple average?", simple)
    message("processing...")
  }
  
  s <- intersect(rownames(countMatrix.center), s)
  message("genes in signature, and also in this dataset: ", length(s))
  
  if (simple){
    sigScore <- colMeans(v[s,])
  } else {
    sigScore <- colMeans(do.call(rbind, lapply(s, function(g) {
      if(verbose) message(".", appendLF = FALSE)
      g.n <- names(sort(abs(countMatrix.mean[g] - countMatrix.mean))[2:(n+1)])
      countMatrix.center[g, ] - colMeans(countMatrix.center[g.n, ])
    })))
  }
  
  if(verbose) message(" done")
  
  metaData[[geneSetName]] <- sigScore
  if (splitBy == "Quartile") { n = 4 }
  if (splitBy == "Median") { n = 2 }
  newColName <- paste0(geneSetName, "Cluster")
  metaData <- metaData %>%
    mutate(paste0(geneSetName, "_Cluster") = ntile(geneSetName, n))
  return(metaData)
}
