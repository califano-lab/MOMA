#' Rank Normalize P Values
#' 
#' This function uses an empirical cumulative distribution function to rank 
#' normalize all the p-values from a particular test
#' @param interacterations.df dataframe off all the interactions
#' @return dataframe will adjusted p values
#' @keywords internal
rankNormalize <- function(interactions.df) {
  # split into two tables, one with event/regulatur/type information, 
  # other with p-values to be ranked/merged. transform these then remerge
  vars <- c("regulator", "event", "type", "sign")
  info.df <- interactions.df[,vars]
  values.df <- interactions.df[,!colnames(interactions.df) %in% vars]
  values.df <- as.matrix(apply(values.df, 2, cdfPval))
  merged.df <- cbind(info.df, values.df)
  merged.df
}

#' Empirical Cumulative Distribution Helper function
#' 
#' This is a wrapper function for the ecdf function in order to get a vector of 
#' the resulting values
#' @importFrom stats ecdf
#' @importFrom tidyr replace_na
#' @param vals original list of p-values to convert
#' @param na.value value to substitute for NAs if present. Default is NA
#' @return vector of adjusted p-values
#' @keywords internal
cdfPval <- function(vals, na.value = NA) {
  # could tweak this...
  vals.adj <- tidyr::replace_na(vals, na.value)
  fn <- stats::ecdf(tidyr::replace_na(vals, 1))
  res <- fn(vals.adj)
  res
}

#' Two step rank integration
#' 
#' This function does the two step integration of the available genomic information
#' First Regulator-Event information is integrated across different tests
#' Then all information for a particular Regulator is integrated to result in a single value
#' that represents the accumulated information for that Regulator
#' @param interactions.df dataframe of all interactions
#' @param na.value value to substitute for NAs if present. Default is NA
#' @return list object with two components: Updated interactions dataframe and a consolidated regulator ranks dataframe
#' @keywords internal
twoStepRankIntegration <- function(interactions.df, na.value) {
  vars <- c("regulator", "event", "type", "sign")
  values.df <- interactions.df[,!colnames(interactions.df) %in% vars]
  
  # Modify first step if NAs are present because otherwise rows 
  # with only 1 value don't need to be integrated and will error 
  message("Doing regulator-event integrations...")
  message("Using ", na.value, " for non existing interaction values.")
  if(!is.na(na.value)) {
    values.df[is.na(values.df)] <- na.value
    event.int.p <- apply(values.df, 1, function(x){
      poolr::fisher(x)$p
    })
  } else {
    event.int.p <- apply(values.df, 1, function(x) {
      if(sum(!is.na(x)) > 1) {
        poolr::fisher(x[!is.na(x)])$p
      } else {
        x[!is.na(x)]
      }
    })
  }
  
  interactions.df$int.p <- event.int.p
  
  # second ranking step
  message("Doing integration of all events for each regulator...")
  ranks.df <- interactions.df %>% 
    dplyr::group_by(.data$regulator) %>%
    dplyr::summarize(int.mr.p = poolr::fisher(.data$int.p)$p)
  
  ranks.df$int.mr.p <- stats::p.adjust(ranks.df$int.mr.p, method = "fdr")
  ranks.df$ecdf.p <- cdfPval(ranks.df$int.mr.p)
  
  return(list(interactions.df = interactions.df, 
              ranks.df = ranks.df))
  
}


#' Over-representation Test for genomic events
#' 
#' Do a proportion test for representation of a particular genomic event in a cluster
#' 
#' @param mat binary events matrix
#' @param inCluster.samples names of samples in the cluster to be tested
#' @param min.events.per.cluster minimum number of times that event must occur in that cluster to be countable
#' @return names of genes that had marginal over-representation 
#' @keywords internal
overrep.analysis <- function(mat, inCluster.samples, min.events.per.cluster) {
  
  # proportion test for over/under representation
  prop.pvals <- apply(mat, 1, function(row, inCluster.samples) {
    row <- na.omit(row)
    if (all(row==0)) { return (1) }
    
    iC.vals <- na.omit(row[inCluster.samples])
    oC.vals <- na.omit(row[setdiff(names(row), inCluster.samples)])
    
    pval <- NULL
    # consider only 0,1 values	
    # look for over-representation
    tab <- rbind( 
      c(length(which(iC.vals>0)), length(which(iC.vals==0)) ),
      c(length(which(oC.vals>0)), length(which(oC.vals==0)) ) 
    )
    rownames(tab) <- c("in-cluster", "out-cluster")
    
    # set threshold for minimum number of times the event 
    # must occur in that cluster
    if (tab[1,1] < min.events.per.cluster) { 
      pval <- 1 
    } else {
      pval <- fisher.test(tab, alternative = "greater")$p.value    
    }
    #pval <- prop.test(tab, alternative='greater')$p.value
    
    pval
    
  }, inCluster.samples=inCluster.samples)
  # anything marginally over represented
  names(prop.pvals[which(prop.pvals < 0.5)])
}
