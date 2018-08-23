#!/usr/bin/env	Rscript

source('moma/R/moma.r')
source('moma/R/tcga.surv.r')

# load data/analyses
load('test_input/viper.rda')
load('test_input/gbm-cindy.rda')
load('test_input/gbm-rawcnv.rda')
load('test_input/gbm-rawsnp.rda')
load('test_input/preppi.pvals.ENTREZ.rda')
load('test_input/gbm-fCNV.rda')

fusions <- as.matrix(read.table("test_input/gbm.txt", header=T, check.names=F, row.names=1))

gene.loc.mapping <- read.table("test_input/metadata/name.loc.map.txt", header=T)

#  CINDY and PREPPI used here for ranking
pathways = list()
pathways[['cindy']] = cindy
pathways[['preppi']] = pval.map

momaObj <- moma.constructor(vipermat, rawsnp, rawcnv, fusions, pathways, 
	gene.blacklist='test_input/mutSig_blacklist.entrezID.txt', 
	output.folder='test_output/',
	gene.loc.mapping=gene.loc.mapping)

momaObj$runDIGGIT(fCNV=fCNV)
momaObj$makeInteractions(cindy.only=FALSE)
momaObj$Rank(use.cindy=TRUE)

# cluster
clustering.solutions <- momaObj$Cluster()
# use clinical survival data to select the clustering solution / break ties in analytical solution. 
# set the clustering variable to MOMA object
clinical <- get.clin('test_input/gdac.broadinstitute.org_GBM.Merge_Clinical.Level_1.2016012800.0.0/GBM.clin.merged.txt')
res <- get.best.clustering.supervised(cluster.sweep=clustering.solutions, clinical=clinical, tissue=tissue, progression.free.surv=TRUE)
momaObj$sample.clustering <- res$clustering
# clustering -> pass to genomic coverage. 
#  
save(momaObj, file='momaObj.rda')

momaObj$saturationPlots()

res <- sapply(1:5, function(x) {
	print (paste("Running iteration ", x))

	set.seed(x)
	viper.submat <- vipermat[,sample(colnames(vipermat), length(colnames(vipermat))*.8)]

	momaObj <- moma.constructor(viper.submat, rawsnp, rawcnv, fusions, pathways, 
		gene.blacklist='test_input/mutSig_blacklist.entrezID.txt', 
		output.folder='test_output/',
		gene.loc.mapping=gene.loc.mapping)
	momaObj$runDIGGIT(fCNV=fCNV)
	momaObj$makeInteractions(cindy.only=FALSE)
	momaObj$Rank(use.cindy=TRUE)

	# return rank orders
	momaObj$ranks$integrated
})
