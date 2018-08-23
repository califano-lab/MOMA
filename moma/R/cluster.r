		


cluster.sweep <- function() { 

	dist.obj <- as.dist(viperSimilarity(w.vipermat, ws=2, method="two.sided"))
	search.results <- clusterRange(dist.obj, range=as.numeric(c(opt$start, opt$end)), step=1, cores=2, method='pam')

}

