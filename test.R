


ScoreSignatureAndSplit <- function(countMatrix.center, countMatrix.mean, s, geneSetName, metaData,
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
  
  sigScore <- sigScore %>% as.data.frame()
  colnames(sigScore)[1] <- geneSetName
  
  if (splitBy == "Quartile") { n = 4 }
  if (splitBy == "Median") { n = 2 }
  
  scoreSplit <- sigScore %>% mutate(clusters = ntile(sigScore[[geneSetName]], n))
  colnames(scoreSplit)[2] <- paste0(geneSetName, "_Cluster")
  
  metaData <- cbind(metaData, scoreSplit)
  
  return(metaData)
}
