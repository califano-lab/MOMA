
##
## Interface definitions section:
##
library(R6)

source("moma/R/cnv.functions.r")
source("moma/R/diggit.r")
library(qvalue)

momaRunner <- R6Class("momaRunner", 
	public = list(
		viper="matrix",
		mut="matrix",
		cnv="matrix",
		fusions="matrix",
		pathways="list",
		gene.blacklist="character",
		output.folder="character",
		cytoband.map="character",
		nes="list", # result field
		interactions="list", # result field
		ranks="list",
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
		
			# need to save the exact hypotheses we're testing, based on MIN.EVENTS
			amps.hypotheses <- rownames(amps[apply(amps,1,sum,na.rm=TRUE)>=min.events,])
			amps.hypotheses <- amps.hypotheses[which(!(amps.hypotheses %in% gene.blacklist))]
		
			dels.hypotheses <- rownames(dels[apply(dels,1,sum,na.rm=TRUE)>=min.events,])
			dels.hypotheses <- dels.hypotheses[which(!(dels.hypotheses %in% gene.blacklist))]
		
			mut.hypotheses <- rownames(somut[apply(somut,1,sum,na.rm=TRUE)>=min.events,])
			mut.hypotheses <- mut.hypotheses[which(!(mut.hypotheses %in% gene.blacklist))]
		
			fus.hypotheses <- rownames(fusions[apply(fusions,1,sum,na.rm=TRUE)>=min.events,])
		
			# Save aREA results
			# Save aREA results
			print ("Writing hypotheses...")
			write.table(amps.hypotheses,file=paste0(output.folder, "/hypotheses.amps.txt"), quote=F, sep="\t")
			write.table(dels.hypotheses,file=paste0(output.folder, "/hypotheses.dels.txt"), quote=F, sep="\t")
			write.table(mut.hypotheses,file=paste0(output.folder, "/hypotheses.muts.txt"), quote=F, sep="\t")
			write.table(fus.hypotheses,file=paste0(output.folder, "/hypotheses.fusions.txt"), quote=F, sep="\t")
		
			# do aREA association
			nes.amps <- associate.events(viper, amps, min.events=min.events, blacklist=gene.blacklist)
			nes.dels <- associate.events(viper, dels, min.events=min.events, blacklist=gene.blacklist)
			nes.muts <- associate.events(viper, somut, min.events=min.events, blacklist=gene.blacklist)
			nes.fusions <- associate.events(viper, fusions, min.events=min.events, blacklist=gene.blacklist)
			
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
			sig.tfs <- viper.getSigTFS(viper)
			null.TFs <- setdiff(rownames(viper), sig.tfs)
			print (paste("Found ", length(sig.tfs), " significant TFs from VIPER scores ... "))
			print (paste("Building background model from ", length(null.TFs), " nulls..."))
	
			local.interactions = list()
			for (type in genomic.event.types) {
				print (paste("Performing background correction for ", type, " NES scores..."))
		 		corrected.scores <- get.diggit.empiricalQvalues(viper, nes[[type]], null.TFs)
				print ("Generating final interactions...")
				local.interactions[[type]] <- sig.interactors.DIGGIT(corrected.scores, nes[[type]], pathways[['cindy']], cindy.only=FALSE)
			}
	
			interactions <<- local.interactions
		},
		Rank = function() {

			# ranks from DIGGIT scores
			z.amp <- cnvScoreStouffer(cytoband.map, interactions[['amp']])
			z.del <- cnvScoreStouffer(cytoband.map, interactions[['del']])
			z.fusions <- stouffer.integrate.diggit(interactions[['fus']])
			z.mut <- stouffer.integrate.diggit(interactions[['mut']])


		}
)


# constructor
moma.constructor <- function(viper, mut, cnv, fusions, pathways, gene.blacklist=NULL, output.folder=NULL, cytoband.mapping=NULL) {

	viper <- samplename.filter(viper)
	mut <- samplename.filter(mut)
	cnv <- samplename.filter(cnv)

	# parse blacklist file, save gene ids
	gene.blacklist = as.character(read.table(gene.blacklist, header=F)[,1])

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
	print (paste("Number of samples in VIPER + Mutation data:", length(nVM))) 
	if (length(nVM) < 2) {
		print ("Error: VIPER and mutation matrix samples don't match!")
		q();
	}
	nVM <- intersect(colnames(viper), colnames(cnv))
	print (paste("Number of samples in VIPER + CNV data:", length(nVM)))
	if (length(nVM) < 2) {
		print ("Error: VIPER and CNV matrix samples don't match!")
		q();
	}

	samples.common <- intersect(intersect(colnames(viper), colnames(mut)), colnames(cnv))
	print (paste(" Common samples for analysis: ", length(samples.common)))
	viper <- viper[,samples.common]
	cnv <- cnv[,samples.common]
	mut <- mut[,samples.common]

	# check TF rows in the indexes	
	for (pathway in names(pathways)) {
		print (paste("Checking labels on pathway", pathway))
		I <- intersect(rownames(viper), names(pathways[[pathway]]))
		print (paste("Found labels for ", length(I), " TFs in VIPER matrix"))
	}
 
	obj <- momaRunner$new(viper=viper, mut=mut, cnv=cnv, fusions=fusions, pathways=pathways, gene.blacklist=gene.blacklist, output.folder=output.folder, cytoband.map=cytoband.mapping)
}

