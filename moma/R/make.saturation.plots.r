#' Helper function to get subtype specific events
#' @param momaObj : object generated from main moma functions
#' @param saturation.data : genomic saturation object from moma. List indexed by cluster then sample then regulator with the number of events associated with each additional regulator
#' @param sample.clustering
#' @return a table that has counts of how many times a particular event happens in a cluster
#' @export
get.subtype.event.tables <- function(saturation.data, sample.clustering) {
  
  subtype.coverage <- list()
  coverage.allSubtypes <- saturation.data
  clustering = sample.clustering
  
  # aggregate statistics across all samples, but use different subtype specific solutions
  for (cluster.id in unique(clustering)) {
    # dataframe with sample, type, gene names for each 
    # amp/del events can be locations. 
    coverage <- coverage.allSubtypes[[cluster.id]]
    sample.set <- names(coverage)
    
    # get number of regulators decided for this checkpoint depending on coverage threshold 
    if (isTRUE(checkpoint.spec)) {
      mr.cutoff <- length(momaObj$checkpoints.clusterSpecific[[cluster.id]])
    } else {
      mr.cutoff <- length(momaObj$checkpoints[[cluster.id]])
    }
    
    subtype.coverage[[cluster.id]] <- make.coverage.df(coverage, cutoff=mr.cutoff)
    
  }
  names(subtype.coverage) <- unique(clustering)
  
  # make summary table that counts the number of events in the group
  subtype.tables <- lapply(names(subtype.coverage), function(cluster.id) {
    table(apply(subtype.coverage[[cluster.id]], 1, function(x) paste0(x[1], ':', x[2])))
  })
  subtype.tables
}
