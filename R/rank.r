#' Use Stouffer's method to combine z-scores of DIGGIT interactions for each cMR protein.
#' 
#' This function combines only positively associated DIGGIT scores by default to create a culmulative DIGGIT score for each cMR.
#' 
#' @import stats
#' @param interactions A list indexed by TF, includes z-scores or p-values for each interacting event
#' @param from.p Integrate p-values or z-scores (default z-scores; from.p = FALSE)
#' @param pos.nes.only Use only positive NES scores to rank proteins (default TRUE)
#' @return A list indexed by TF, a stouffer integrated z-score
stoufferIntegrateDiggit <- function(interactions, from.p = FALSE, pos.nes.only = TRUE) {
  
  ## Integrate p-values for each TF
  diggit.integrated.z <- unlist(lapply(interactions, function(scores) {
    # stouffer's method
    if (from.p) {
      scores <- -qnorm(scores)
    }
    
    if (length(scores) == 1) {
      if (pos.nes.only) {
        if (as.numeric(scores) > 0) {
          return(as.numeric(scores))
        } else {
          return(0)
        }
      }
    }
    
    integrated.z <- 0
    if (pos.nes.only) {
      scores <- scores[which(as.numeric(scores) > 0)]
      integrated.z <- sum(scores)/sqrt(length(scores))
      if (is.nan(integrated.z)) {
        integrated.z <- 0
      }
    } else {
      integrated.z <- sum(abs(scores))/sqrt(length(scores))
    }
    integrated.z
  }))
  names(diggit.integrated.z) <- names(interactions)
  diggit.integrated.z
}


#' @title dispatch method for either CNV location corrected or SNV
#' @param interactions List of MR - Genomic Event interactions, inferred by DIGGIT
#' @param cytoband.map Data.frame mapping Entrez.IDs to cytoband locations
#' @return Z-scores for each MR
stoufferIntegrate <- function(interactions, cytoband.map = NULL) {
  z <- NULL
  if (!is.null(cytoband.map)) {
    # need to create a vector with gene
    map.vec <- cytoband.map$Cytoband
    names(map.vec) <- cytoband.map$Entrez.IDs
    z <- cnvScoreStouffer(map.vec, interactions)
  } else {
    z <- stoufferIntegrateDiggit(interactions)
  }
  z
}


#' Integrate CNV scores 
#' 
#' @param mapping a named vector of genomic locations/cytoband IDs. 
#' names are the gene names for each--i.e. a many to one mapping from HUGO or entrez IDs to cytoband location
#' @param diggit.interactions list indexed by MR/TF name in Entrez Space each points to a named vector of NES / z-scores associated with entrez IDs for each interacting event. 
#' @param cytoband Boolean to use cytoband locations for computing final integrated score
#' @param from.p Boolean, set TRUE if diggit.interaction values are p-values instead of z-scores
#' @param pos.nes.only Boolean, only consider positive DIGGIT association scores when ranking candidate MRs (default=TRUE)
#' @return A vector of z-scores, named by the Master Regulators in 'diggit.interactions' 
cnvScoreStouffer <- function(mapping, diggit.interactions, cytoband = TRUE, from.p = FALSE, pos.nes.only = TRUE) {
  
  mapped.diggit <- mapScoresCnvBand(mapping, diggit.interactions, from.p = FALSE, pos.nes.only = TRUE)
  
  # apply over each TF/MR: sum the named vector of scores
  integrated.z.scores <- unlist(lapply(mapped.diggit, function(scores) {
    
    # now integrate scores with Stouffer's method
    integrated.z <- sum(scores)/sqrt(length(scores))
    if (is.null(integrated.z)) {
      integrated.z <- 0
    }
    if (is.nan(integrated.z)) {
      integrated.z <- 0
    }
    integrated.z
  }))
  
  names(integrated.z.scores) <- names(diggit.interactions)
  integrated.z.scores
}

#' Map scores to cytoband location 
#' 
#' @param mapping a named vector of genomic locations/cytoband IDs. 
#' names are the gene names for each--i.e. a many to one mapping from HUGO or entrez IDs to
#' cytoband location
#' @param diggit.interactions list indexed by MR/TF name in Entrez Space
#' @param from.p DIGGIT interactions are in p-value format instead of z-score (default=FALSE)
#' @param pos.nes.only Only consider positive associations with NES scores (default=TRUE)
#' each points to a named vector of NES / z-scores associated with entrez IDs for each interacting event.
#' @return A list of input scores, now named by cytoband location
mapScoresCnvBand <- function(mapping, diggit.interactions, from.p = FALSE, pos.nes.only = TRUE) {
  
  # apply over each TF/MR:
  mapped.scores <- lapply(diggit.interactions, function(x) {
    
    # print (paste('Integrating: ', length(x), ' interactions ')) check for zero interactions case
    if (length(x) == 0) {
      return(0)
    }
    
    scores <- NULL
    # corner case of 1 interaction only...
    if (length(x) == 1) {
      score <- 0
      if (from.p) {
        score <- qnorm(min(x), lower.tail = FALSE)
      } else {
        if (pos.nes.only) {
          # Only include positive scores
          score <- x
          if (score < 0) {
            score <- 0
          }
        } else {
          # Only include positive scores
          score <- x
        }
      }
      names(score) <- names(x)
      scores <- score
    } else {
      
      # gather either cytoband locations or Locus names are the genes in ENTREZ space
      
      locations <- na.omit(mapping[as.character(names(x))])
      
      if (length(locations) == 0) {
        print(paste("Warning: didn't find any genomic locations for", names(x)))
        return(0)
      }
      
      # find the highest MR-gene score in this location and use that as a summary
      scores <- unlist(lapply(unique(locations), function(loc) {
        # for each location, get entrez gene IDs
        genes.this.loc <- as.character(names(locations[which(locations == loc)]))
        # determine if we're using z-scores or p-values here
        scores.thisLoc <- x[genes.this.loc]
        if (length(scores.thisLoc) == 0) {
          return(NA)
        }
        score <- NULL
        if (from.p) {
          score <- qnorm(min(scores.thisLoc), lower.tail = FALSE)
        } else {
          if (pos.nes.only) {
            # Only include positive scores
            score <- as.numeric(max(scores.thisLoc))
            if (score < 0) {
              score <- 0
            }
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


#' @title Combine DIGGIT inferences with pathway knowledge 
#' @param diggit.int List of interactions between MRs - Genomic events, inferred by DIGGIT
#' @param pathway - a list indexed by TF/MR entrez ID, contains the named vector of p-values for interactions 
#' @param pos.nes.only Only use positive associations between MR activity and presence of events (default = True)
#' @param cores Number of cores to use if parallel is selected
#' @return numeric vector, zscores for each TF/MR
pathwayDiggitIntersect <- function(diggit.int, pathway, pos.nes.only = TRUE, cores = 1) {
  
  
  pathway.pvals <- parallel::mclapply(names(diggit.int), function(tf) {
    
    partners.pvals <- pathway[[as.character(tf)]]
    
    # 
    if (pos.nes.only) {
      I <- diggit.int[[as.character(tf)]]
      I <- I[which(I > 0)]
      tf.diggit.interactors <- unique(names(I))
    } else {
      tf.diggit.interactors <- unique(names(diggit.int[[as.character(tf)]]))
    }
    
    # Find partners PrePPI P-values
    pvals <- partners.pvals[which(names(partners.pvals) %in% tf.diggit.interactors)]
    
    if (length(pvals) == 0) {
      return(1)
    }
    if (length(pvals) == 1) {
      return(as.numeric(pvals))
    }
    # print (paste0('found this many relevant interactions: ', length(pvals))) print (pvals)
    pvals
  }, mc.cores = cores)
  names(pathway.pvals) <- names(diggit.int)
  
  ## Compute both Stouffer integrated z-scores and Fisher integrated p-values
  integrated.z <- unlist(lapply(pathway.pvals, function(x) {
    
    x[which(x == 1)] <- 0.5
    iz <- sum(abs(qnorm(x)))/sqrt(length(x))
    iz
  }))
  names(integrated.z) <- names(pathway.pvals)
  
  integrated.p <- unlist(lapply(pathway.pvals, function(x) {
    
    if (length(x) == 1) {
      return(x)
    }
    
    integrated.p <- -pchisq(-2 * sum(log(x)), 2 * length(x), log.p = TRUE)
    integrated.p
  }))
  names(integrated.p) <- names(diggit.int)
  
  # return the z-scores only
  integrated.z
  
}


#' Implements the conditional Bayes model to combine VIPER scores with diggit and pathway scores
#'
#' @param viper.scores numeric Vector 
#' @param pathway.scores List , double indexed by each pathway dataset, then with type char. Each points to a numeric score vectors in [0,R+] for each
#' @param diggit.scores List indexed by type char, with numeric score vectors in [0,R+] for each
#' @return a named vector of empirical p-values for each protein/candidate Master Regulator
conditionalModel <- function(viper.scores, diggit.scores, pathway.scores) {
  
  integrated.pvals <- unlist(lapply(names(viper.scores), function(VIP) {
    # VIP = Viper Inferred Protein
    
    # VIPER scores are the independent model component
    viper.p <- empiricalP(VIP, viper.scores)
    
    # for each type, produce p-values for DIGGIT as well as each dependent pathway or pathway based algorithm (including CINDy)
    all.pvals <- c(viper.p)
    for (type in names(diggit.scores)) {
      
      # DIGGIT is conditioned on VIPER scores
      diggit.p <- conditionalP(VIP, viper.scores, diggit.scores[[type]])
      # Pathway scores are conditioned on DIGGIT
      pathway.pvals <- c()
      for (pathway in names(pathway.scores)) {
        if (length(pathway.scores[[pathway]]) == 0) {
          next
        } else if (!(type %in% names(pathway.scores[[pathway]]))) {
          stop(paste("Error: can't find entry for ", type, " in pathway ranks ", pathway))
        }
        pathway.p <- conditionalP(VIP, diggit.scores[[type]], pathway.scores[[pathway]][[type]])
        pathway.pvals <- c(pathway.p, pathway.pvals)
      }
      all.pvals <- unlist(c(all.pvals, diggit.p, pathway.pvals))
    }
    
    x <- na.omit(all.pvals)
    # Fisher's method
    integrated.p <- 1 - pchisq(-2 * sum(log(x)), 2 * length(x), log.p = FALSE)
    integrated.p
  }))
  names(integrated.pvals) <- names(viper.scores)
  integrated.pvals
}


#' Get the empirical p-value from a distribution (vector)
#' 
#' @param gene.name Character
#' @param x named Vector of scores for the distribution
#' @return a numeric p-value between 0 and 1
empiricalP <- function(gene.name, x) {
  
  # unranked genes in either distribution should get no score from this
  if (!(as.character(gene.name) %in% names(x))) {
    return(1)
  }
  
  n.gte <- length(which(rank(x) >= rank(x)[as.character(gene.name)]))
  p.one.tail <- n.gte/(length(x) + 1)
  p.one.tail
}


#' Get the conditional p-value of a gene
#'  
#' @param gene.name Character
#' @param condition.on named Vector of scores for the distribution we are conditioning ON
#' @param x named Vector of scores for the dependent distribution
#' @return a numeric p-value between 0 and 1
conditionalP <- function(gene.name, condition.on, x) {
  
  # null ranks for all: return NA
  if (all(x == 0)) {
    return(NA)
  }
  
  # unranked genes in either distribution should get no score from this
  if (!(as.character(gene.name) %in% names(condition.on))) {
    return(1)
  }
  if (!(as.character(gene.name) %in% names(x))) {
    return(1)
  }
  
  # these names are GTE the gene of interest
  gte <- names(which(rank(condition.on) >= rank(condition.on)[as.character(gene.name)]))
  
  conditional.dist <- x[gte]
  
  n.gte <- length(which(rank(conditional.dist) >= rank(conditional.dist)[as.character(gene.name)]))
  # empirical p-value: note the GTE means we count the value itself
  p.one.tail <- n.gte/(length(conditional.dist) + 1)
  p.one.tail
  
}



