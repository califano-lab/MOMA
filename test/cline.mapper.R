
# STES test dataset

library(moma)
source('moma/R/models.match.r')

load('data/STES_vipermat.rda')
load('data/ccle.large_intensine.rda')

patient.vipermat <- vipermat
model.vipermat <- vpmat

# test using a single checkpoint and single cluster that includes all patients
STAD.checkpoint.entrezID <- as.character(read.table('data/stad.cMRs.txt', header=T)[,1])
checkpoints <- list()
checkpoints[[1]] <- STAD.checkpoint.entrezID
clusters <- list()
clusters[[1]] <- colnames(patient.vipermat)

matchObj <- viperMatch.constructor(patient.viper=vipermat, 
			model.viper=vpmat, 
			checkpoints=checkpoints,
			clusters=clusters)

matchObj$checkpoint.enrichments()
matchObj$similarity()
