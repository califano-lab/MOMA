#' Generate clusters from viper values 
#'
#' @import clusterpam
#' @import parallel
#' @export
cluster.sweep <- function() { 

	dist.obj <- as.dist(viper::viperSimilarity(w.vipermat, ws=2, method="two.sided"))
	search.results <- clusterpam::clusterRange(dist.obj, range=as.numeric(c(opt$start, opt$end)), step=1, cores=2, method='pam')

}

