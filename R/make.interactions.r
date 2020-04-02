#' Utility function
#' 
#' @param vipermat - matrix of VIPER scores with columns as samples, rows as protein names
#' @param fdr.thresh - BH-FDR threshold (default 0.05 FDR rate)
#' @examples
#' viper.getTFScores(gbm.example$vipermat)
#' @return A vector of normalized z-scores, named by TF id
#' @export
viper.getTFScores <- function(vipermat, fdr.thresh = 0.05) {
  
  # for each gene, count the number samples with scores for each, and weight by that
  w.counts <- apply(vipermat, 1, function(x) {
    data.counts <- length(which(!is.na(x)))
    data.counts
  })
  w.counts <- w.counts/ncol(vipermat)
  
  vipermat[is.na(vipermat)] <- 0
  
  # normalize element scores to sum to 1 (optional - use weighted element scores based on silhouette)
  element.scores <- rep(1, ncol(vipermat))
  element.scores <- element.scores/sum(element.scores)
  
  # mean weighted VIPER score across samples
  w.means <- apply(vipermat, 1, function(x) {
    res <- sum(x * element.scores)
    res
  })
  # weight by the counts for each
  w.means <- w.means * w.counts
  names(w.means) <- rownames(vipermat)
  
  # only look at those with positive (high) score w.means <- sort(w.means[which(w.means > 0)], decreasing=TRUE)
  
  zscores <- w.means
  zscores
}


#' Calculate p-values from pseudo zscores / VIPER aREA scores, threshold
#' 
#' @param zscores Vector of normally distributed z-scores representing protein activities. 
#' @param fdr.thresh Threshold for false discovery rate, default is 0.05
#' @return Get the names of proteins with significant z-scores, after multi-hypothesis correction
viper.getSigTFS <- function(zscores, fdr.thresh = 0.05) {
  
  # calculate pseudo-pvalues and look at just significant pvals/scores
  pvals <- -pnorm(abs(zscores), log.p = TRUE) * 2
  pvals[which(pvals > 1)] <- 1
  # correct unless option is NULL
  sig.idx <- which(p.adjust(pvals, method = "BH") < fdr.thresh)
  pvals <- pvals[sig.idx]
  
  names(pvals)
}


#' Compute the empirical q-values of each genomic-event/VIPER gene pair 
#' 
#' Use against the background distribution of associations with a given set of 'null' VIPER genes (i.e. low activity TFs) 
#' 
#' @param vipermat viper inferences matrix, samples are columns, rows are TF entrez gene IDs 
#' @param nes scores for each mutation (rows) against each TF (columns) 
#' @param null.TFs low-importance TFs used to calculate null distributions
#' @param alternative Alternative defaults to 'both' : significant p-values can come from both sides of the null distribution 
#' @return A named list of qvalues for each TF/cMR protein. Each entry contains a vector of q-values for all associated events; names are gene ids 
get.diggit.empiricalQvalues <- function(vipermat, nes, null.TFs, alternative = "both") {
  
  # subset NES to Viper Proteins in the vipermat only
  nes <- nes[, as.character(rownames(vipermat))]
  
  nes.em.qvals <- apply(nes, 1, function(x, alternative) {
    
    null.VEC <- x[as.character(null.TFs)]
    null.VEC <- null.VEC[which(!is.na(null.VEC))]
    # get empirical q-values for both upper and lower tails of NES / DIGGIT statistics
    qvals <- get.empirical.qvals(x, null.VEC, alternative)
    qvals
  }, alternative = alternative)
  
  names(nes.em.qvals) <- rownames(nes)
  nes.em.qvals
}


#' Get.empirical.qvals
#' 
#' @param test.statistics P-values generated from the test comparisons
#' @param null.statistics P-values generated under the null (permutation) model
#' @param alternative Optional : 1 or 2 tails used to generate the p-value (default='both')
#' @return A list with both the qvalues and empirical p-values from the supplied test and null stats
get.empirical.qvals <- function(test.statistics, null.statistics, alternative = "both") {
  
  # calculate the upper and lower tail
  if (alternative == "both") {
    
    test.statistics <- sort(abs(test.statistics), decreasing = TRUE)
    null.statistics <- abs(null.statistics)
    
    em.pvals <- qvalue::empPvals(test.statistics, null.statistics)
    qvals <- rep(1, length(em.pvals))
    tryCatch({
      qvals <- qvalue::qvalue(em.pvals)$qvalue
    }, error = function(e) {
      # if pi0, the estimated proportion of true null hypothesis <= 0, it might fail: in that case set to zero and return p-values anyways
      qvals <- rep(1, length(em.pvals))
    })
    names(qvals) <- names(test.statistics)
    names(em.pvals) <- names(test.statistics)
    return(list(qvals = qvals, pvals = em.pvals))
  } else {
    stop(paste(" alternative ", alternative, " not implemented yet!"))
  }
}


#' Filter interactions from NES (DIGGIT) scores and corresponding background-corrected scores. 
#' 
#' Use this version in the Bayes model to rank TFs
#' 
#' @import stats
#' @param corrected.scores A list indexed by the genomic event/gene with corresponding pvals and qvals for
#' each TF
#' @param nes.scores Matrix with tfs as columns, rows are genomic events
#' @param cindy CINDy algorithm output matrix
#' @param p.thresh P-value threshold (default=0.05)
#' @param cindy.only Consider only CINDy validated interactions (default=TRUE)
#' @return a list (indexed by VIPER protein) of significant genomic interactions 
#' and associated pvals over the background (null TF) model, and NES scores
sig.interactors.DIGGIT <- function(corrected.scores, nes.scores, cindy, p.thresh = 0.05, cindy.only = TRUE) {
  
  pvals.matrix <- get.pvals.matrix(corrected.scores)
  
  # input validation
  if (!is.numeric(p.thresh)) {
    print("Error: invalid value supplied for p-value threshold!")
    q()
  }
  
  ## Apply joint NES + p-value over background (null TF) threshold over each Viper Protein return the raw NES scores only for those significant over the
  ## background and including CINDy (if applicable)
  viper.interactors <- lapply(colnames(pvals.matrix), function(viperProt) {
    
    # find the over-null-TF-background scores with an significant, uncorrected p-value print (paste('Processing : ', viperProt))
    pvals <- as.numeric(pvals.matrix[, as.character(viperProt)])
    nes.vec <- as.numeric(nes.scores[, as.character(viperProt)])
    # print (nes.vec) subset to significant p-values
    row.idx <- which(pvals < p.thresh)
    pvals <- pvals[row.idx]
    names(pvals) <- rownames(pvals.matrix)[row.idx]
    
    # if (!all(names(pvals) == names(nes.vec))) { print ('Error: data not aligned for aREA / aREA corrected p-values') }
    
    # subset the NES vector for this TF, threshold again on NES scores as a sanity check on the basic enrichment (i.e. remove those with high
    # over-background scores simply because the background is de-enriched)
    nes.vec <- nes.scores[, which(colnames(nes.scores) == as.character(viperProt))]
    nes.vec <- nes.vec[which(2 * (1 - pnorm(abs(nes.vec))) < p.thresh)]
    
    fusion.index <- unlist(lapply(names(nes.vec), function(x) if (length(strsplit(x, "_")[[1]]) > 1) 
      TRUE else FALSE))
    # subset to CINDY validated upstream regulators
    if (cindy.only && is.null(cindy[[viperProt]])) {
      return(c())
    } else if (cindy.only) {
      
      # lookup takes about a second...
      cindy.thisTF <- cindy[[viperProt]]
      
      upstream.cindy.modulators <- names(cindy.thisTF)
      cindy.nes.vec <- na.omit(nes.vec[match(upstream.cindy.modulators, names(nes.vec))])
      if (length(cindy.nes.vec) == 0) {
        return(c())
      }
      
      # add fusions in, they don't have CINDy scores
      fus.vec <- nes.vec[fusion.index]
      entrez.vec <- cindy.nes.vec[which(!is.na(cindy.nes.vec))]
      
      # Independence of these: multiply through
      cindy.entrez.pvals <- cindy.thisTF[(names(cindy.thisTF) %in% names(entrez.vec))]
      entrez.pvals <- 2 * (1 - pnorm(abs(entrez.vec)))
      entrez.pvals <- entrez.pvals * cindy.entrez.pvals
      
      entrez.vec.corrected <- (1 - qnorm(entrez.pvals)) * sign(entrez.vec)
      nes.vec <- c(entrez.vec.corrected, fus.vec)
    }
    
    # keep interactions significant above the null model and with a significant raw aREA score
    nes.vec <- nes.vec[intersect(names(nes.vec), names(pvals))]
    nes.vec <- nes.vec[which(!is.na(nes.vec))]
    
    # print (paste('num interactions:', length(nes.vec)))
    nes.vec <- sort(nes.vec, decreasing = TRUE)
    nes.vec
  })
  names(viper.interactors) <- colnames(pvals.matrix)
  viper.interactors
}

