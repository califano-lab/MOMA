#!/usr/bin/env	Rscript

library(moma)
library(parallel)
library(clusterpam)
library(magrittr)

# load data/analyses
load('test/gbm/viper.rda')
load('test/gbm/gbm-cindy.rda')
load('test/gbm/gbm-rawcnv.rda')
load('test/gbm/gbm-rawsnp.rda')
load('test/preppi.pvals.ENTREZ.rda')
load('test/gbm/gbm-fCNV.rda')

fusions <- as.matrix(read.table("test/gbm/gbm-fusions_PRADA-calls.txt", header=T, check.names=F, row.names=1))

gene.loc.mapping <- read.table("test/name.loc.map.txt", header=T)

#  CINDY and PREPPI used here for ranking
pathways = list()
pathways[['cindy']] = cindy
pathways[['preppi']] = pval.map

momaObj <- moma.constructor(vipermat, rawsnp, rawcnv, fusions, pathways, 
	gene.blacklist='test/mutSig_blacklist.entrezID.txt', 
	output.folder='gbm-test/',
	gene.loc.mapping=gene.loc.mapping)

momaObj$runDIGGIT(fCNV=fCNV)
momaObj$makeInteractions(cindy.only=FALSE)
momaObj$Rank(use.cindy=TRUE)

# cluster
clustering.solutions <- momaObj$Cluster()
# use clinical survival data to select the clustering solution / break ties in analytical solution. 
# set the clustering variable to MOMA object
clinical <- get.clin('test/gbm/GBM.clin.merged.txt')
res <- get.best.clustering.supervised(cluster.sweep=clustering.solutions, clinical=clinical, tissue=tissue, progression.free.surv=TRUE)
momaObj$sample.clustering <- res$clustering
# clustering -> pass to genomic coverage. 
#  

#save(momaObj, file='momaObj.rda')

momaObj$saturationPlots()

res <- sapply(1:5, function(x) {
	print (paste("Running iteration ", x))

	set.seed(x)
	viper.submat <- vipermat[,sample(colnames(vipermat), length(colnames(vipermat))*.8)]

	momaObj <- moma.constructor(viper.submat, rawsnp, rawcnv, fusions, pathways, 
		gene.blacklist='test/mutSig_blacklist.entrezID.txt', 
		output.folder='gbm-test/',
		gene.loc.mapping=gene.loc.mapping)
	momaObj$runDIGGIT(fCNV=fCNV)
	momaObj$makeInteractions(cindy.only=FALSE)
	momaObj$Rank(use.cindy=TRUE)

	# return rank orders
	momaObj$ranks$integrated
})

save(momaObj, file='momaObj.rda')

