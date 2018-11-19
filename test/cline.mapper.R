
# STES test dataset

library(moma)

load('data/ccle.CNS.rda')
load('momaObj.gbm.rda')

patient.vipermat <- momaObj$viper
model.vipermat <- vpmat

checkpoints <- momaObj$checkpoints
clusters <- momaObj$sample.clustering

matchObj <- moma::viperMatch.constructor(patient.viper=patient.vipermat, 
			model.viper=model.vipermat, 
			checkpoints=checkpoints,
			clusters=lapply(unique(clusters), function(x) names(clusters[clusters==x])))

matchObj$checkpoint.enrichments()
matchObj$similarity()
