
library(reshape)
library(methods)

#' 
#' Map scores to cytoband location 
#' 
#' @param mapping: a named vector of genomic locations/cytoband IDs. 
#' names are the gene names for each--i.e. a many to one mapping from HUGO or entrez IDs to
#' cytoband location
#' @param diggit.interactions: list indexed by MR/TF name in Entrez Space
#' each points to a named vector of NES / z-scores associated with entrez IDs for each interacting event.
#' @export
mapScores.cnvBand <- function(mapping, diggit.interactions, from.p=FALSE, pos.nes.only=TRUE) {

	# apply over each TF/MR:
	mapped.scores <- lapply(diggit.interactions, function(x) {

		#print (paste("Integrating: ", length(x), " interactions "))
		# check for zero interactions case
		if (length(x) == 0) {
			return (0)
		}	

		scores <- NULL
		# corner case of 1 interaction only...
		if (length(x) == 1) {
			score <- 0
			if (from.p) {
				score <- qnorm(min(x), lower.tail=F)
			} else {
				if (pos.nes.only) {
					# Only include positive scores
					score <- x
					if (score < 0) { score <- 0 }
				} else {
					# Only include positive scores
					score <- x
				}
			}
			names(score) <- names(x)
			scores <- score
		} else {

			# gather either cytoband locations or Locus
			# names are the genes in ENTREZ space
			# 
		
			locations <- na.omit(mapping[as.character(names(x))])

			if (length(locations) == 0) { 
				print (paste("Warning: didn't find any genomic locations for", names(x)))
				return (0)
			}
		
			# find the highest MR-gene score in this location and use that as a summary
			scores <- unlist(lapply(unique(locations), function(loc) {
				# for each location, get entrez gene IDs
				genes.this.loc <- as.character(names(locations[which(locations == loc)]))
				# determine if we're using z-scores or p-values here
				scores.thisLoc <- x[genes.this.loc]
				if (length(scores.thisLoc) == 0) { return (NA) }
				score <- NULL
				if (from.p) {
					score <- qnorm(min(scores.thisLoc), lower.tail=F)
				} else {
					if (pos.nes.only) {
						# Only include positive scores
						score <- as.numeric(max(scores.thisLoc))
						if (score < 0) { score <- 0 }
					} else {
						# take the max abs value of the score
						score <- scores.thisLoc[which.max(abs(scores.thisLoc))]
					}
				}
				score
			}))
			names(scores) <- unique(locations)
		}

		# add fusion events if necessary
		scores
	})
	names(mapped.scores) <- names(diggit.interactions)
	mapped.scores
}


#' 
#' Integreate CNV scores 
#' 
#' @param mapping: a named vector of genomic locations/cytoband IDs. 
#' names are the gene names for each--i.e. a many to one mapping from HUGO or entrez IDs to
#' cytoband location
#' @param diggit.interactions: list indexed by MR/TF name in Entrez Space
#' 	each points to a named vector of NES / z-scores associated with entrez IDs for each interacting event. 
#' @export
cnvScoreStouffer <- function(mapping, diggit.interactions, cytoband=TRUE, from.p=FALSE, pos.nes.only=TRUE) {

	mapped.diggit <- mapScores.cnvBand(mapping, diggit.interactions, from.p=FALSE, pos.nes.only=TRUE)

	# apply over each TF/MR: sum the named vector of scores 
	integrated.z.scores <- unlist(lapply(mapped.diggit, function(scores) {

		# now integrate scores with Stouffer's method
		integrated.z <- sum(scores)/sqrt(length(scores))
                if (is.null(integrated.z)) { integrated.z <- 0 }
                if (is.nan(integrated.z)) { integrated.z <- 0 }
		integrated.z
	}))

	names(integrated.z.scores) <- names(diggit.interactions)
	integrated.z.scores
}

