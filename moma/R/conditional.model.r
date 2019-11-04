

#'
#' @title Get the empirical p-value from a distribution (vector)
#' @param gene.name Character
#' @param x named Vector of scores for the distribution
#' @return a numeric p-value between 0 and 1
empirical.p <- function(gene.name, x) {
    
    # unranked genes in either distribution should get no score from this
    if (!(as.character(gene.name) %in% names(x))) {
        return(1)
    }
    
    n.gte <- length(which(rank(x) >= rank(x)[as.character(gene.name)]))
    p.one.tail <- n.gte/(length(x) + 1)
    p.one.tail
}

#'
#' @title Get the conditional p-value of a gene 
#' @param gene.name Character
#' @param condition.on named Vector of scores for the distribution we are conditioning ON
#' @param x named Vector of scores for the dependent distribution
#' @return a numeric p-value between 0 and 1
conditional.p <- function(gene.name, condition.on, x) {
    
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

#' Implements the conditional Bayes model to combine VIPER scores with diggit and pathway scores
#'
#' @param viper.scores numeric Vector 
#' @param pathway.scores List , double indexed by each pathway dataset, then with type char. Each points to a numeric score vectors in [0,R+] for each
#' @param diggit.scores List indexed by type char, with numeric score vectors in [0,R+] for each
#' @return a named vector of empirical p-values for each protein/candidate Master Regulator
conditional.model <- function(viper.scores, diggit.scores, pathway.scores) {
    
    integrated.pvals <- unlist(lapply(names(viper.scores), function(VIP) {
        # VIP = Viper Inferred Protein
        
        # VIPER scores are the independent model component
        viper.p <- empirical.p(VIP, viper.scores)
        
        # for each type, produce p-values for DIGGIT as well as each dependent pathway or pathway based algorithm (including CINDy)
        all.pvals <- c(viper.p)
        for (type in names(diggit.scores)) {
            
            # DIGGIT is conditioned on VIPER scores
            diggit.p <- conditional.p(VIP, viper.scores, diggit.scores[[type]])
            # Pathway scores are conditioned on DIGGIT
            pathway.pvals <- c()
            for (pathway in names(pathway.scores)) {
                if (length(pathway.scores[[pathway]]) == 0) {
                  next
                } else if (!(type %in% names(pathway.scores[[pathway]]))) {
                  stop(paste("Error: can't find entry for ", type, " in pathway ranks ", pathway))
                }
                pathway.p <- conditional.p(VIP, diggit.scores[[type]], pathway.scores[[pathway]][[type]])
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


