#!/usr/bin/env	Rscript

library(moma)
#library(parallel)
#library(magrittr)

# load data/analyses
load('test/gbm/gbm.data.rda')

#  CINDY and PREPPI used here for ranking
pathways = list()
pathways[['cindy']] = cindy
pathways[['preppi']] = pval.map

mutsig.blacklist <- as.character(read.table('test/mutSig_blacklist.entrezID.txt', header=F)[,1])

momaObj <- moma.constructor(vipermat, rawsnp, rawcnv, fusions, pathways, 
	gene.blacklist=mutsig.blacklist,
	output.folder='gbm-test/',
	gene.loc.mapping=gene.loc.mapping)

momaObj$runDIGGIT(fCNV=fCNV)
momaObj$makeInteractions(cindy.only=FALSE)
momaObj$Rank(use.cindy=TRUE)

# cluster
clustering.solutions <- momaObj$Cluster()
# Pick the k=5 cluster and save this to the moma Object
momaObj$sample.clustering <- clustering.solutions[[4]]$clustering
# generate genomic saturation statistics and plot data
momaObj$saturationPlots()

cluster1.checkpoint <- momaObj$checkpoints[[1]]

