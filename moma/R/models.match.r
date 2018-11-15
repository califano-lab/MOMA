
##
## Functions to match cell line and other model systems or cross-patient
## groups.
##

#' @import viper
#' @import tidyverse
#' @import RColorBrewer
#' @export
model_viperMatcher <- setRefClass("model_viperMatcher", fields=
	list(patient.viper="matrix",
	model.viper="matrix",
	checkpoints="list",
	clusters="list",
	vpsim="matrix",
	CP.enrichments="matrix"),
	methods = list(
		#' 
		#' @param checkpoints : a list with each index containing a checkpoint definition (entrezIDS) to match 
		#' against a set of model system VIPER scores
		#' @param model.vipermat : viper matrix for the model system being tested. Columns are samples, rows proteins
		#' 
		checkpoint.enrichments = function() {
			checkpoint.enrichments <- moma::aREA.regulon_enrich(checkpoints, t(model.viper))
			rownames(checkpoint.enrichments) <- names(checkpoints)
			CP.enrichments <<- checkpoint.enrichments
		},
		similarity = function() {
			joint.matrix <- as.matrix(cbind(patient.viper, model.viper))
			vpsim <<- as.matrix(viperSimilarity(joint.matrix, ws=c(4, 2), method='two.sided'))
		}
	)
)

viperMatch.constructor <- function(patient.vipermat, model.vipermat, checkpoints, clusters) {

	patient.vipermat <- patient.vipermat[intersect(rownames(patient.vipermat), rownames(model.vipermat)),]
	model.vipermat <- model.vipermat[intersect(rownames(patient.vipermat), rownames(model.vipermat)),]

	for (cp in names(checkpoints)) {
		checkpoint <- checkpoints[[cp]]
		print (paste("Warning, matching checkpoint with ", 
			length(which(checkpoint %in% rownames(model.vipermat))), " of ", 
			length(checkpoint), " checkpoint proteins in the cell line matrix.."))
		for (missing.prot in setdiff(checkpoint, rownames(model.vipermat))) {
			print (paste("Missing protein: ", map.entrez(missing.prot)))
		}
	}

	obj <- model_viperMatcher$new(patient.viper=patient.vipermat, model.viper=model.vipermat, checkpoints=checkpoints, clusters=clusters)
	obj
}

