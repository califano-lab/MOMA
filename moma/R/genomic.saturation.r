suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))

#' Get coverage of interactions
#' 
#' @param momaObj A numeric vector with cluster membership, names are samples
#' @param viper.samples Calculate the genomic coverage only for these sample
#' @param cMR.ranking A vector entrez IDs, in order
#' @param topN Compute coverage for only the top -N- Master Regulators
#' @param mutation.filter Retain only mutation events in this (positive) list
#' @export
get.coverage <- function(momaObj, cMR.ranking, viper.samples, topN=100, mutation.filter=NULL) {

	if (class(momaObj) != 'momaRunner') {
		stop("Error: must have instantiated momaRunner class object passed!")
	}

	# select considered cMRs
	print (paste("Top : ", topN, " regulators"))
	selected.tfs <- cMR.ranking[1:topN]
	if (length(selected.tfs)==0) {
		print ("Error: no TFs selected!")
		q();
	}
	
	# interactions
	# makeInteractions = function(genomic.event.types=c("amp", "del", "mut", "fus"), cindy.only=TRUE) { 
	# interactions: 
	# 	mut : named by Entrez ID
	interaction.map <- valid.diggit.interactions(momaObj$interactions, momaObj$gene.loc.mapping, selected.tfs)
	
	# another assert statment: make sure we have non-zero interactions for each
	sapply(names(interaction.map), function(key) {
		if (sum(sapply(interaction.map[[key]], function(x) length(x))) < 1) {
			print(paste("Warning: didn't find any positive DIGGIT associations for data type ", key))
			print(paste("(in subtype)"))
		}
	})

	oc <- sample.overlap(momaObj, viper.samples, selected.tfs, interaction.map, mutation.filter=mutation.filter)
	# count mutations/amps/dels covered at this point. Aggregate stats

	oc
}

#' The core function to compute which sample-specific alterations overlap with genomic events that are explained 
#' via DIGGIT. 
#' 
#' @param momaObj Object reference of momaRunner class
#' @param viper.samples Sample vector to restrict sample-specific analysis to
#' @param selected.tfs Transcription factors being analyzed
#' @param interaction.map List object of events 'covered' by the supplied interactions of type mut/amp/del/fus
#' @param mutation.filter A vector of whitelisted mutation events, in entrez gene IDs
#' @param idx.range Number of tfs to check for genomic saturation calculation, default is 1253
#' @param cnv.threshold Numeric absolute value to threshold SNP6 and/or GISTIC or other CNV scores at. Above that absolute value is considered a positive event. 
#' @param verbose Output status during the run (default=FALSE)
#' @export
#' 
sample.overlap <- function(momaObj, viper.samples, selected.tfs, interaction.map, cnv.threshold=0.5, mutation.filter=NULL, verbose=FALSE, idx.range=NULL) {

	if (is.null(momaObj$hypotheses)) {
		stop('Error: no hypothesis set for momaRunner class object!!')
	}

	map <- momaObj$gene.loc.mapping

	mut.HYP.filter <- momaObj$hypotheses[['mut']]
	amp.HYP.filter <- na.omit(unique(sapply(as.character(momaObj$hypotheses[['amp']]), function(x) {
		res <- NA
		x <- as.character(x)
		if (x %in% map$Entrez.IDs) {
			# get cytoband location for this gene if it's in the map
			res <- unique(map[which(map$Entrez.IDs==x),"Cytoband"])
		} 
		as.character(res)
	})))
	del.HYP.filter <- na.omit(unique(sapply(as.character(momaObj$hypotheses[['del']]), function(x) {
		res <- NA
		x <- as.character(x)
		if (x %in% map$Entrez.IDs) {
			# get cytoband location for this gene if it's in the map
			res <- unique(map[which(map$Entrez.IDs==x),"Cytoband"])
		} 
		as.character(res)
	})))
	fus.HYP.filter <- momaObj$hypotheses[['fus']]

	if (length(mut.HYP.filter)==0) {
		stop('Error: null hypotheses for mut!')
	} else if (length(amp.HYP.filter)==0) {
		stop('Error: null hypotheses for amp!')
	} else if (length(del.HYP.filter)==0) {
		stop('Error: null hypotheses for del!')
	} else if (length(fus.HYP.filter)==0) {
		warning("Zero fusion hypotheses")
	}
	# consider only samples triple intersection of samples with CNV, mutation and RNA (i.e. VIPER inferences) data
	all.samples.genomics <- intersect(colnames(momaObj$mut), colnames(momaObj$cnv))
	all.sample.names <- intersect(all.samples.genomics, viper.samples)
	coverage <- lapply(all.sample.names, function(sample) {

		print (paste0("Computing coverage for sample ",  sample))

		# find active and inactive proteins in this sample
		viper.pvals <- (1-pnorm(abs(momaObj$viper[,sample])))
		active.proteins <- rownames(momaObj$viper)[intersect(which(viper.pvals < 0.05), which(momaObj$viper[,sample] > 0))]
		if (length(active.proteins)==0) {
			warning(paste0("No active proteins found for sample: ", sample))
		}

		# Collect mutation events in this sample's row:
		mut.events <- as.character(rownames(momaObj$mut[which(momaObj$mut[,sample] > 0),]))
		mut.events <- mut.events[which(mut.events %in% momaObj$hypotheses[['mut']])]

		# entrez ids to cytoband
		# get genes
		del.events <- as.character(rownames(momaObj$cnv[which(momaObj$cnv[,sample] < -cnv.threshold),]))
		del.events.entrez <- del.events[which(del.events %in% momaObj$hypotheses[['del']])]
		# map to genomic locations
		del.events.cytoband <- unique(map[which(map$Entrez.IDs %in% del.events.entrez),"Cytoband"])

		# get genes
		amp.events <- as.character(rownames(momaObj$cnv[which(momaObj$cnv[,sample] > cnv.threshold),]))
		amp.events.entrez <- amp.events[which(amp.events %in% momaObj$hypotheses[['amp']])]
		# map to genomic locations
		amp.events.cytoband <- unique(map[which(map$Entrez.IDs %in% amp.events.entrez),"Cytoband"])

		# not all samples will have fusions: include if possible
		fus.events <- NULL
		if (!is.null(fus.HYP.filter)) {
			if (sample %in% colnames(momaObj$fus)) {
				fus.events <- rownames(momaObj$fus[which(momaObj$fus[,sample] > 0),])
			} else {
				warning(paste0("Sample not recorded in fusions matrix:", sample))
			}
		}
		fus.events <- fus.events[which(fus.events %in% momaObj$hypotheses[['fus']])]

		validated.del.locations <- del.events.cytoband
		validated.amp.locations <- amp.events.cytoband
		validated.fusion.events <- fus.events
		validated.mut.events <- mut.events

		if (!is.null(mutation.filter)) {
			print ("Using mutation filter:")
			prev.count <- length(validated.mut.events)
			validated.mut.events <- intersect(mutation.filter, validated.mut.events)
			removed = prev.count - length(validated.mut.events)
			print (paste("Filtered out: ", removed, " mutations using filter..."))
		}

		if (verbose) {
			print (validated.fusion.events)
			print ('AMPS:')
			print (validated.del.locations)
			print ('DELS:')
			print (validated.amp.locations)
			print ('Muts:')
			print (validated.mut.events)
		}

		# for each K in 1:N, compute the coverage of the top K events
		# and return the stats
		sample.cover <- lapply(1:length(selected.tfs), function(x) c())
		# do a semi-complete range for speed		
		# and if supplied already with a range, then use that
		if (is.null(idx.range)) {
			if (length(selected.tfs) > 100) {
				idx.range <- 1:50
				idx.range <- c(idx.range, 26:50*2)
				idx.range <- c(idx.range, 11:30*10)
				idx.range <- c(idx.range, 13:25*25)
				idx.range <- c(idx.range, 7:12*100)
				idx.range <- c(idx.range, 1253)
			} else {
				idx.range = 1:length(selected.tfs)
			}
		} 

		# store the set of active MRs, if it's identical to the last iteration, update and skip re-computation
		active.mrs.length.LAST <- 0
		last.idx <- 1
		for (i in idx.range) { 
			# use these for coverage
			top.N.tfs <- selected.tfs[1:i]
			active.mrs <- intersect(active.proteins, top.N.tfs)

			# if no extra active mrs this round, copy previous stats and continue
			if ((length(active.mrs)==active.mrs.length.LAST) & (i>1)) {
				sample.cover[[i]] <- sample.cover[[last.idx]]
				next
			} else {
				active.mrs.length.LAST <- length(active.mrs)
				last.idx <- i
			}

			#print ("Active MRs:")
			#print (active.mrs)
			covered.mut <- unique(unlist(lapply(active.mrs, function(mr) { covered <- intersect(interaction.map$mut[[mr]], validated.mut.events); covered })))
			covered.del <- unique(unlist(lapply(active.mrs, function(mr) { covered <- intersect(interaction.map$del[[mr]], validated.del.locations); covered })))
			covered.amp <- unique(unlist(lapply(active.mrs, function(mr) { covered <- intersect(interaction.map$amp[[mr]], validated.amp.locations); covered })))
			covered.fus <- unique(unlist(lapply(active.mrs, function(mr) { covered <- intersect(interaction.map$fus[[mr]], validated.fusion.events); covered })))
			amp.frac <- ifelse(length(validated.amp.locations)==0, NA, length(covered.amp)/length(validated.amp.locations))
			del.frac <- ifelse(length(validated.del.locations)==0, NA, length(covered.del)/length(validated.del.locations))
			mut.frac <- ifelse(length(validated.mut.events)==0, NA, length(covered.mut)/length(validated.mut.events))
			fus.frac <- ifelse(length(validated.fusion.events)==0, NA, length(covered.fus)/length(validated.fusion.events))

			# compute total frac
			total.frac <- NA 
			total.num.events <- length(c(validated.mut.events, validated.amp.locations, validated.del.locations, validated.fusion.events))

			total.frac <- ifelse(length(total.num.events)==0, NA, (length(covered.amp)+length(covered.del)+length(covered.mut)+length(covered.fus))/total.num.events)
			sample.cover[[i]] <- list(amp=covered.amp, 
						del=covered.del, 
						mut=covered.mut, 
						fus=covered.fus, 
						amp.frac=amp.frac, 
						del.frac=del.frac, 
						mut.frac=mut.frac, 
						fus.frac=fus.frac, 
				total.frac=total.frac)
		}
		sample.cover
	})
	names(coverage) <- all.sample.names
	coverage
}

#' Return a set of events 'covered' by specified cMR-event interactions 
#' @param interactions List indexed by amp/mut/del/fus - from cMRs to interacting events
#' @param selected.tfs For each event type list, search within only these cMRS
#' @param gene.loc.mapping Data.frame mapping entrezIDs to cytoband locations
#' @return a list of events 'covered' by the supplied interactions of type mut/amp/del/fus
#' @export
valid.diggit.interactions <- function(interactions, gene.loc.mapping, selected.tfs) {
	
	if (length(selected.tfs)==0) {
		stop("No TFs input to diggit function")
	}
	selected.tfs <- as.character(selected.tfs)

	mut.tfs <- selected.tfs[which(selected.tfs %in% names(interactions[['mut']]))]
	amp.tfs <- selected.tfs[which(selected.tfs %in% names(interactions[['amp']]))]
	del.tfs <- selected.tfs[which(selected.tfs %in% names(interactions[['del']]))]
	fus.tfs <- selected.tfs[which(selected.tfs %in% names(interactions[['fus']]))]

	if (length(mut.tfs)==0 || length(amp.tfs)==0 || length(del.tfs)==0) {
		stop("No valid TFs in the interactions supplied!")
	} else if (length(fus.tfs)==0) {
		warning("No valid fusion interactions...")
	}

	mut.I <- subset.list.interactions(interactions[['mut']], mut.tfs)
	del.I <- subset.list.interactions(interactions[['del']], del.tfs)
	amp.I <- subset.list.interactions(interactions[['amp']], amp.tfs)
	
	# only subset if fusions exist
	# fusions IDs are unique , nothing else necessary
	if (length(fus.tfs) >= 1){
	  fus.I <- subset.list.interactions(interactions[['fus']], fus.tfs)
	  covered.fusions <- fus.I
	}
	

	# Add cnv events to mutation coverage, as either copy number variation is 
	# valid evidence for explaining a patient's mutation  
	covered.mutations <- merge.lists(mut.I, del.I)
	covered.mutations <- merge.lists(covered.mutations, amp.I)

	# Add mut events to cnv coverage, as it is a valid 
	# type of evidence for explaining a patient's CNV change  
	covered.amps <- merge.lists(amp.I, mut.I)
	covered.dels <- merge.lists(del.I, mut.I)

	

	# create a new mapping from TF in Entrez -> event location
	covered.amps.LOC <- lapply(names(covered.amps), function(x, I) {
		geneNames <- I[[as.character(x)]]
		band.names <- unique(as.character(gene.loc.mapping[which(gene.loc.mapping$Entrez.IDs %in% geneNames),"Cytoband"]))
		if (length(band.names)==0) {
			print(paste("Warning: could not map entrez IDs to Cytoband for IDS, skipping..."))
			print(geneNames)
			band.names <- NA
		}
		band.names
	}, I=covered.amps)
	names(covered.amps.LOC) <- names(covered.amps)
	covered.amps.LOC <- covered.amps.LOC[!is.na(covered.amps.LOC)]

	if (sum(sapply(covered.amps.LOC, function(x) length(x)))==0) {
		print("Error: something went wrong when mapping amplification Entrez.IDs to Cytoband IDs. Quitting...")
		quit(status=1)
	}

	# create a new mapping from TF in Entrez -> event location
	covered.dels.LOC <- lapply(names(covered.dels), function(x, I) {
		geneNames <- I[[as.character(x)]]
		band.names <- unique(as.character(gene.loc.mapping[which(gene.loc.mapping$Entrez.IDs %in% geneNames),"Cytoband"]))
		if (length(band.names)==0) {
			print(paste("Warning: could not map entrez IDs to Cytoband for IDS, skipping..."))
			print(geneNames)
			band.names <- NA
		}
		band.names
	}, I=covered.dels)
	names(covered.dels.LOC) <- names(covered.dels)
	covered.dels.LOC <- covered.dels.LOC[!is.na(covered.dels.LOC)]

	if (sum(sapply(covered.amps.LOC, function(x) length(x)))==0) {
		print("Error: something went wrong when mapping deletion Entrez.IDs to Cytoband IDs. Quitting...")
		quit(status=1)
	}
	
	# don't incorporate fusions into final list object unless they exist 
	if (length(fus.tfs) >= 1) {
	  return (list(mut=covered.mutations, amp=covered.amps.LOC, del=covered.dels.LOC, fus=covered.fusions))
	} else {
	  return (list(mut=covered.mutations, amp=covered.amps.LOC, del=covered.dels.LOC)) }
}

#' @title Helper function: subset a list to the set of keys supplied return the names of interactions with positive values, in a list structure
#' @param int.l List of interactions, at each index this is a numeric named vector
#' @param keys Keys used to reduce interactions
subset.list.interactions <- function(int.l, keys) {

	filtered.I <- lapply(keys, function(key, interactions) {
		I <- interactions[[as.character(key)]]
		# named vector
		event.names <- names(I[which(I > 0)])
		event.names	
	}, interactions=int.l)
	names(filtered.I) <- keys

	filtered.I
}

merge.lists <- function(l1, l2) {

	merged <- list()
	inter <- intersect(names(l1), names(l2))
	# combine joint entries
	for (key in inter) {
		merged[[key]] <- union(l1[[key]], l2[[key]])
	}	
	# add entries unique to l1
	for (key in setdiff(names(l1), names(l2))) {
		merged[[key]] <- l1[[key]]
	}
	# add entries unique to l2
	for (key in setdiff(names(l2), names(l1))) {
		merged[[key]] <- l2[[key]]
	}
	return (merged)
}
	
