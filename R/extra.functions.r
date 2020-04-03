#' Compute aREA enrichment for the proteins in a given regulon, against vipermat scores supplied
#' (used for model matching functionality)
#' 
#' @param regulon ARACNE regulon object
#' @param vipermat A VIPER network of inferred activity scores with columns as patient samples, and rows as proteins
#' @return A matrix of enrichment scores with rows as event/gene names and columns as VIPER protein names
aREA.regulon_enrich <- function(regulon, vipermat) {
  
  regulon <- lapply(regulon, function(x) as.character(x))
  # Calculate raw enrichment scores: each mutation against each TF columns are TFs, rownames are mutations
  es <- rea(t(vipermat), regulon)
  # Analytical correction
  dnull <- reaNULL(regulon)
  
  # Calculate pvalue of ES
  pval <- t(vapply(seq_len(length(dnull)), function(i, es, dnull) {
    dnull[[i]](es[i, ])$p.value
  }, es = es$groups, dnull = dnull, FUN.VALUE = numeric(ncol(es$groups))))
  
  # Convert the pvalues into Normalized Enrichment Scores
  nes <- qnorm(pval/2, lower.tail = FALSE) * es$ss
  # print (dim(nes)) print (length(regulon))
  rownames(nes) <- names(regulon)
  colnames(nes) <- rownames(vipermat)
  dimnames(pval) <- dimnames(nes)
  nes[is.na(nes)] <- 0
  # columns are TFs, rows are genomic events
  nes
}


#' Simple one-tail rank based enrichment analysis sREA
#' (for cluster analysis)
#' 
#' This function performs simple 1-tail rank based enrichment analysis
#' 
#' @param signatures Numeric matrix of signatures
#' @param groups List containing the groups as vectors of sample names
#' @return Matrix of Normalized Enrichment Zcores
sREA <- function(signatures, groups) {
  if (is.null(nrow(signatures))) 
    signatures <- matrix(signatures, length(signatures), 1, dimnames = list(names(signatures), "sample1"))
  # ranked signatures matrix: samples are rows. Rank
  sig <- qnorm(apply(signatures, 2, rank)/(nrow(signatures) + 1))
  gr <- sapply(groups, function(x, samp) {
    samp %in% x
  }, samp = rownames(sig))
  gr <- t(gr)
  # non negative counts
  nn <- rowSums(gr)
  # fractions of prev computation rows are groups
  gr <- gr/nn
  es <- gr %*% sig
  return(es * sqrt(nn))
}

#' Utility function 
#' 
#' @param corrected.scores - corrected p-values processed by 'qvals' package
#' @return A matrix of p-values for scores between genes/events (rows) and TFs (columns)
get.pvals.matrix <- function(corrected.scores) {
  # order of VIPER proteins/TFs
  tf.names.order <- names(corrected.scores[[1]]$qvals)
  pvals.matrix <- matrix(unlist(lapply(corrected.scores, function(x) {
    pvals <- x$pvals[tf.names.order]
    pvals
  })), byrow = TRUE, ncol = length(tf.names.order))
  
  colnames(pvals.matrix) <- tf.names.order
  rownames(pvals.matrix) <- names(corrected.scores)
  pvals.matrix
}


#' Retain TCGA sample ids without the final letter designation ('A/B/C') 
#' 
#' @param input Matrix of expression or protein activity scores. Columns are sample names, rows are genes.
#' Input can also just be an input vector of sample names.
#' @examples 
#' sample.names <- c("TCGA-14-1825-01A", "TCGA-76-4931-01B", "TCGA-06-5418-01A")
#' samplename.filter(sample.names)
#' @return An identical matrix with new (shorter) column names, or a vector with the shortened names. 
#' @export
samplename.filter <- function(input) {
  # filter down to sample Id without the 'A/B/C sample class'.
  if(is.matrix(input)) {
    sample.ids <- vapply(colnames(input), function(x) substr(x, 1, 15), FUN.VALUE = character(1))
    colnames(input) <- sample.ids
  } else if (is.vector(input)) {
    input <- substr(input, 1, 15)
  }
  
  input
}




