
##
## Interface definitions section:
##

#source("R/cnv.functions.r")
#source("R/diggit.r")
#source("R/conditional.model.r")
#source("R/genomic.saturation.r")

#' @import qvalue
#' @import parallel
#' @import atools
#' @import clusterpam
#' @import reshape2
#' @import tidyverse
#' @import MKmisc
#' @import survival
#' @import RColorBrewer
#' @export
momaRunner <- setRefClass("momaRunner", fields=
	list(viper="matrix",
	mut="matrix",
	cnv="matrix",
	fusions="matrix",
	pathways="list",
	gene.blacklist="character",
	output.folder="character",
	gene.loc.mapping="data.frame",
	nes="list", # result field
	interactions="list", # result field
	clustering.results="list",
	ranks="list",
	hypotheses="list",
	genomic.saturation="list",
	sample.clustering="numeric"), # numbers are cluster assignments, names are sample ids matching other data
	methods = list(
		runDIGGIT = function(fCNV=NULL, cnvthr=0.5, min.events=4) {
	
			cnv.local <- NULL
			if (is.null(fCNV)) {
				print ("Warning: no fCNV supplied, using no CNV filter!")
				cnv.local <- cnv
			} else {
				cnv.local <- cnv[intersect(fCNV,rownames(cnv)),]
			}
		
			somut <- mut
		
			# Define amplifications and deletions
			amps <- dels <- cnv.local
			amps[which(amps < cnvthr)] <- 0
			amps[which(amps >= cnvthr)] <- 1
			dels[which(dels > -cnvthr)] <- 0
			dels[which(dels <= -cnvthr)] <- 1
		
			# save the exact hypotheses (genes) we're testing, based on MIN.EVENTS
			amps.hypotheses <- rownames(amps[apply(amps,1,sum,na.rm=TRUE)>=min.events,])
			amps.hypotheses <- amps.hypotheses[which(!(amps.hypotheses %in% gene.blacklist))]
		
			dels.hypotheses <- rownames(dels[apply(dels,1,sum,na.rm=TRUE)>=min.events,])
			dels.hypotheses <- dels.hypotheses[which(!(dels.hypotheses %in% gene.blacklist))]
		
			mut.hypotheses <- rownames(somut[apply(somut,1,sum,na.rm=TRUE)>=min.events,])
			mut.hypotheses <- mut.hypotheses[which(!(mut.hypotheses %in% gene.blacklist))]
		
			# Save aREA results
			# Save aREA results
			print ("Writing hypotheses...")
			dir.create(output.folder, showWarnings=FALSE)
			write.table(amps.hypotheses,file=paste0(output.folder, "/hypotheses.amps.txt"), quote=F, sep="\t")
			write.table(dels.hypotheses,file=paste0(output.folder, "/hypotheses.dels.txt"), quote=F, sep="\t")
			write.table(mut.hypotheses,file=paste0(output.folder, "/hypotheses.muts.txt"), quote=F, sep="\t")
			hypotheses <<- list('mut' = mut.hypotheses, 'del' = dels.hypotheses, 'amp' = amps.hypotheses)	

			# do aREA association
			nes.amps <- associate.events(viper, amps, min.events=min.events, blacklist=gene.blacklist)
			nes.dels <- associate.events(viper, dels, min.events=min.events, blacklist=gene.blacklist)
			nes.muts <- associate.events(viper, somut, min.events=min.events, blacklist=gene.blacklist)

			nes.fusions <- NULL
			if (!is.null(fusions)) {
				fus.hypotheses <- rownames(fusions[apply(fusions,1,sum,na.rm=TRUE)>=min.events,])
				write.table(fus.hypotheses,file=paste0(output.folder, "/hypotheses.fusions.txt"), quote=F, sep="\t")
				nes.fusions <- associate.events(viper, fusions, min.events=min.events, blacklist=gene.blacklist)
			}
			
			# Save aREA results
			save(nes.amps, nes.dels, nes.muts, nes.fusions, file=paste0(output.folder, "/aREA.rda"))
			# store in the object list
			nes <<- list(amp=nes.amps, del=nes.dels, mut=nes.muts, fus=nes.fusions)
		},
		makeInteractions = function(genomic.event.types=c("amp", "del", "mut", "fus"), cindy.only=TRUE) { 

			# NULL TFs :
			# from the VIPER matrix, calculate p-values of the absolute mean NES score for each. 
			# (2-tailed test, -pnorm*2). 
			# BH-FDR < 0.05 are sig. Take everything else as the background model. 
			ranks[['viper']] <<- viper.getTFScores(viper)
			sig.tfs <- viper.getSigTFS(ranks[['viper']])
			null.TFs <- setdiff(rownames(viper), sig.tfs)
			print (paste("Found ", length(sig.tfs), " significant TFs from VIPER scores ... "))
			print (paste("Building background model from ", length(null.TFs), " nulls..."))
	
			local.interactions = list()
			for (type in genomic.event.types) {
				print (paste("Performing background correction for ", type, " NES scores..."))
				nes.thisType <- nes[[type]]
				if (is.null(nes.thisType)) { 
					next 
				}
		 		corrected.scores <- get.diggit.empiricalQvalues(viper, nes.thisType, null.TFs)
				print ("Generating final interactions...")
				local.interactions[[type]] <- sig.interactors.DIGGIT(corrected.scores, nes[[type]], pathways[['cindy']], cindy.only=FALSE)
			}
	
			interactions <<- local.interactions
		},
		Rank = function(use.cindy=FALSE, genomic.event.types=c("amp", "del", "mut", "fus")) {
			# ranks from DIGGIT scores
			integrated.z <- list()
			for (type in genomic.event.types) {
		
				if (is.null(interactions[[type]])) {
					next
				}
				# deletions/amp CNV events need to be corrected for genome location...
				if (type == 'amp' || type == 'del') {
					integrated.z[[type]] <- stouffer.integrate(interactions[[type]], gene.loc.mapping)
				} else {
					integrated.z[[type]] <- stouffer.integrate(interactions[[type]], NULL)
				}
			}
			# generate integrated rankings from additional sources of evidence, including CINDy and
			# pathway databases and/or structural databases like PrePPI
			pathway.z <- list()	
			for (pathway in names(pathways)) {
	
				pathway.z[[pathway]] <- list()
				# optional: use cindy scores to generate a separate ranking	
				if (pathway=='cindy' && !use.cindy) {
					next
				}
				for (type in genomic.event.types) {
					if (is.null(interactions[[type]])) {
						next
					}
					pathway.z[[pathway]][[type]] <- pathway.diggit.intersect(interactions[[type]], pathways[[pathway]], pos.nes.only=TRUE)
				}
			}

			viper.scores <- ranks[['viper']]	
			print ("Integrating all data in Bayesian conditional model...")
			ranks[['integrated']] <<- conditional.model(viper.scores, integrated.z, pathway.z)
		}, 
		Cluster = function() {
			# do weighted pearson correlation, using the ranks as weights
			weights <- log(ranks[['integrated']])^2
		 	weights <- weights[as.character(rownames(viper))]
                	w.vipermat <- weights*viper
                	print ("using pearson correlation with weighted vipermat")
                	dist.obj <- corDist(t(w.vipermat), method="pearson")
			print ("testing clustering options, k = 2..15")
			search.results <- clusterRange(dist.obj, range=as.numeric(c(2,15)), step=1, cores=2, method='pam')
			search.results
		}, 
		saturationPlots = function(clustering.solution=NULL, cov.fraction=0.85) {

			# do saturation analysis with all available data types
			# helper functions in external libs
			# coverage statistics at throughhold sweep

			if (is.null(clustering.solution)) {
				clustering.solution <- sample.clustering
			}
			# get coverage for each subtype
			coverage.subtypes <- list()
			for (clus.id in unique(clustering.solution)) {
				viper.samples <- colnames(viper[,names(clustering.solution[clustering.solution==clus.id])])
				coverage.range <- get.coverage(.self, viper.samples, clustering.solution, topN=100)
				coverage.subtypes[[clus.id]] <- coverage.range
			}
			genomic.saturation <<- coverage.subtypes
		}
	)
)



#
#' MOMA Constructor
#' @param mut : an indicator matrix (0/1) of mutation events with samples as columns and genes as rows
#' @param fusions : an indicator matrix (0/1) of fusion events with samples as columns and genes as rows
#' @param cnv : a matrix of CNV scores (typically SNP6 array scores from TCGA) with samples as columns and genes as rows
#' @param pathways : a named list of lists. Each named list represents interactions between proteins (keys) and their associated partners
#' @param cytoband.mapping : vector of band locations, names are Entrez IDs
#' @return an instance of class momaRunner
#' @export
moma.constructor <- function(viper, mut, cnv, fusions, pathways, gene.blacklist=NULL, output.folder=NULL, gene.loc.mapping=NULL) {
	viper <- samplename.filter(viper)
	mut <- samplename.filter(mut)
	cnv <- samplename.filter(cnv)

	# validate viper matrix
	if (ncol(viper) < 2) {
		print ("Error: fewer than 2 samples in VIPER matrix!")
		q();
	}
	if (nrow(viper) < 2) {
		print ("Error: fewer than 2 rows in VIPER matrix!")
		q();
	}

	# check column overlap
	nVM <- intersect(colnames(viper), colnames(mut))
	mut <- mut[,nVM]
	print (paste("Number of samples in VIPER + Mutation data:", length(nVM))) 
	if (length(nVM) < 2) {
		print ("Error: VIPER and mutation matrix samples don't match!")
		q();
	}
	nVM <- intersect(colnames(viper), colnames(cnv))
	cnv <- cnv[,nVM]
	print (paste("Number of samples in VIPER + CNV data:", length(nVM)))
	if (length(nVM) < 2) {
		print ("Error: VIPER and CNV matrix samples don't match!")
		q();
	}

	if (!is.null(fusions)) {
		nVM <- intersect(colnames(viper), colnames(fusions))
		# redo type conversion to matrix: could have a single row with this datatype
		# so we must re-type it to guard against auto conversion to 'integer'
		fusions <- as.matrix(fusions[,nVM])
		print (paste("Number of samples in VIPER + Fusion data:", length(nVM)))
		if (length(nVM) < 2) {
			print ("Error: VIPER and CNV matrix samples don't match!")
			q();
		}
	}

	if (!is.null(gene.loc.mapping)) {
		# verify the 
		if (!class(gene.loc.mapping)=='data.frame') {
			stop("Error: gene location mapping supplied is not a data.frame!")
		} else if (!("Entrez.IDs" %in% colnames(gene.loc.mapping))) {
			stop("Error: gene location mapping supplied does not have 'Entrez.IDs' attribute!")
		} else if (!("Cytoband" %in% colnames(gene.loc.mapping))) {
			stop("Error: gene location mapping supplied does not have 'Cytoband' attribute!")
		}

	} else {
		print ("Warning: no gene - genomic location mapping provided!")
	}

	
	samples.common <- intersect(intersect(colnames(viper), colnames(mut)), colnames(cnv))
	print (paste(" Common samples with all data analysis: ", length(samples.common)))

	# check TF rows in the indexes	
	for (pathway in names(pathways)) {
		print (paste("Checking labels on pathway", pathway))
		I <- intersect(rownames(viper), names(pathways[[pathway]]))
		if (length(I)==0) {
			stop("No intersection with supplied pathway! Assuming a formatting error and quitting...")
		}
		print (paste("Found labels for ", length(I), " TFs in VIPER matrix"))
	}
 
	obj <- momaRunner$new(viper=viper, mut=mut, cnv=cnv, fusions=fusions, pathways=pathways, 
		gene.blacklist=as.character(gene.blacklist), 
		output.folder=output.folder, 
		gene.loc.mapping = gene.loc.mapping)
	obj
}


#' preppi: 3 columns, partner A, B and the p-value of the interaction
#' @param : pathway - a list indexed by TF/MR entrez ID, contains the 
#' 	named vector of p-values for interactions 
#' @returns numeric vector, zscores for each TF/MR
#' @export
pathway.diggit.intersect <- function(diggit.int, pathway, pos.nes.only=TRUE) {

	pathway.pvals <- mclapply(names(diggit.int), function(tf) {

		partners.pvals <- pathway[[as.character(tf)]]

		# 
		if (pos.nes.only) {
			I <- diggit.int[[as.character(tf)]]
			I <- I[which(I > 0)]
			tf.diggit.interactors <- unique(names(I))
		} else {
			tf.diggit.interactors <- unique(names(diggit.int[[as.character(tf)]]))
		}

		#  Find partners PrePPI P-values
		pvals <- partners.pvals[which(names(partners.pvals) %in% tf.diggit.interactors)]

		if (length(pvals) == 0) {
			return (1)
		}
		if (length(pvals) == 1) {
			return (as.numeric(pvals))
		}
		#print (paste0("found this many relevant interactions: ", length(pvals)))
		#print (pvals)
		pvals
	}, mc.cores=10)
	names(pathway.pvals) <- names(diggit.int)

	##
	## Compute both Stouffer integrated z-scores and Fisher integrated p-values
	##
	integrated.z <- unlist(lapply(pathway.pvals, function(x) {

		x[which(x == 1)] <- 0.5
		iz <- sum(abs(qnorm(x)))/sqrt(length(x))
		iz
	}))
	names(integrated.z) <- names(pathway.pvals)

	integrated.p <- unlist(lapply(pathway.pvals, function(x) {

		if (length(x) == 1) { return (x) }

		integrated.p <- -pchisq( -2*sum(log(x)) , 2*length(x), log.p=TRUE)
		integrated.p
	}))
	names(integrated.p) <- names(diggit.int)

	# return the z-scores only
	integrated.z

	#list(pvals=integrated.p, zscores=integrated.z, pathway.pvals=pathway.pvals)
}

#' dispatch method for either CNV location corrected or SNV
#' @export
stouffer.integrate <- function(interactions, cytoband.map=NULL) {
	z <- NULL
	if (!is.null(cytoband.map)) {
		# need to create a vector with gene
		map.vec <- cytoband.map$Cytoband
		names(map.vec) <- cytoband.map$Entrez.IDs
		z <- cnvScoreStouffer(map.vec, interactions)
	} else {
		z <- stouffer.integrate.diggit(interactions)
	}
	z
}
