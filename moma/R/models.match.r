
##
## Functions to match cell line and other model systems or cross-patient
## groups.
##
library(RColorBrewer)
library(gplots)
fn1 <- colorRampPalette(c("#1E90FF", "#FFFFFF"))
fn2 <- colorRampPalette(c("#FFFFFF", "#FF4500"))

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
	models.match="list",
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
			res <- matrix(viperSimilarity(joint.matrix, ws=c(4, 2), method='two.sided'), ncol=ncol(joint.matrix))
			rownames(res) <- colnames(joint.matrix)
			colnames(res) <- colnames(joint.matrix)
			vpsim <<- res
		}, 
		matchingPatients = function(match.threshold=as.numeric("1e-5")) {

			# for each patient cluster, match 
			# determine sample order
			best.cline.models <- list()
			for (clus.id in 1:length(clusters)) {

				# patient samples
				clus.samples <- clusters[[clus.id]]

				# these cell lines matche the checkpoint
				cp.pvals <- 1-pnorm(matchObj$CP.enrichments[clus.id,])
				cp.matches <- names(cp.pvals[p.adjust(cp.pvals, method='BH') < 0.05])

				if (length(cp.matches)==0) {
					print (paste("Warning, didn't find any matches to the checkpoint for cluster ", clus.id))
					print ("Using viperSimilarity() only!")
					cp.matched <- colnames(model.viper)
				}

				# count the number of patients each model system matches, according to viperSimiliarity()
				# at the specified p-value threshold
				number.patients.matched <- sapply(cp.matches, function(cline) {
					patient.scores <- vpsim[clus.samples, cline]
					pvals <- 1-pnorm(patient.scores)
					length(pvals[which(pvals < match.threshold)])
				})

				ranked.models <- sort(number.patients.matched, dec=T)
				best.cline.models[[clus.id]] <- ranked.models
			}
			models.match <<- best.cline.models
	
		}, 
		#' Plot the viper matrix of the best matching cell lines against the 
		#' 
		heatmap.patientCoverage = function(models.match, output.file, top.clines=3) {

			# features to select
			features <- unlist(checkpoints)

			# scale the VIPER patient matrix for display
			zoom.vipermat <- c()
			for (clus.id in 1:length(clusters)) {
				clus.samples <- clusters[[clus.id]]
				zoom.vipermat <- cbind(zoom.vipermat, patient.viper[features,clus.samples])
			}
			scaled.vipermat <- t(scale(t(zoom.vipermat), center=F))

			data <- c()
			for (clus.id in 1:length(clusters)) {
				clus.samples <- clusters[[clus.id]]
				selected.models <- names(models.match[[clus.id]][1:top.clines])
				data <- cbind(data, scaled.vipermat[features,clus.samples], model.viper[features, selected.models])
			}

			#
			cline.total.counts <- sapply(unique(names(unlist(models.match))), function(cline) {
				sum(sapply(1:6, function(m) {
					r <- 0	
					if (cline %in% names(models.match[[m]])) {
						r <- models.match[[m]][cline]
					} 
					r
				}))
			
			})

			# hugo gene names	
			rownames(data) <- map.entrez(rownames(data))
			
			pdf(output.file, width=16, height=16)
			heatmap.2(data, 
				col=c(fn1(25), fn2(25)), 
				colsep=F, 
				rowsep=F, 
				distfun=corDist, 
				tracecol=NA, 
				Colv=FALSE, 
				Rowv=FALSE, 
				margins=c(10,5), 
				keysize=0.5,
			breaks=seq(-3,3, length.out=51))
			dev.off()

		},
		#' Plot the viper matrix of the best matching cell lines against the MEAN of each cluster
		#' 
		heatmap.patientMean = function(output.file) {

			# features to select
			features <- unlist(checkpoints)
	
			# for each patient cluster, match 
			# determine sample order
			data <- c()
			for (clus.id in 1:length(clusters)) {

				# patient samples
				clus.samples <- clusters[[clus.id]]

				# get mean match scores
				model.scores <- sort(apply(vpsim[clus.samples,colnames(model.viper)], 2, mean), dec=T)

				# significant matches to the checkpoint
				cp.pvals <- 1-pnorm(matchObj$CP.enrichments[clus.id,])
				cp.matches <- names(cp.pvals[p.adjust(cp.pvals, method='BH') < 0.05])
				if (length(cp.matches)==0) {
					print (paste("Warning, didn't find any matches to the checkpoint for cluster ", clus.id))
					print ("Using viperSimilarity() only!")
					selected.models <- names(sort(model.scores, dec=T))[1:3]
				} else {
					# get mean match scores
					model.scores <- sort(apply(vpsim[clus.samples,colnames(model.viper)], 2, mean), dec=T)
					selected.models <- names(sort(model.scores[cp.matches], dec=T)[1:min(3, length(cp.matches))])
					# add the patient clusters and the selected models on the right
				}
				data <- cbind(data, patient.viper[features,clus.samples], model.viper[features, selected.models])
			}

			# hugo gene names	
			rownames(data) <- map.entrez(rownames(data))

			# select best match cell lines


			pdf(output.file, width=16, height=16)
			heatmap.2(data, 
				col=c(fn1(25), fn2(25)), 
				colsep=F, 
				rowsep=F, 
				distfun=corDist, 
				tracecol=NA, 
				Colv=FALSE, 
				Rowv=FALSE, 
				margins=c(10,5), 
				keysize=0.5,
			breaks=seq(-4,4, length.out=51))
			dev.off()

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

