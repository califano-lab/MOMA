#' Use 'aREA' to calculate the enrichment between each genomic event - VIPER inferred protein pair. 
#' 
#' Requires pre-computed VIPER scores and a binary events matrix. Will use only samples in both event and VIPER matrices.
#'  
#' @param vipermat Pre-computed VIPER scores with samples as columns and proteins as rows 
#' @param events.mat Binary 0/1 events matrix with samples as columns and genes or events as rows
#' @param whitelist Only compute associations for events in this list
#' @param min.events Only compute enrichment if the number of samples with these events is GTE to this
#' @param event.type Name of the event type being analyzed
#' @return A matrix of aREA scores, dimensions are nrow(events.mat) x nrow(vipermat) 
associateEvents <- function(vipermat, events.mat, min.events = NA, whitelist = NA, event.type = c("Amplifications", "Deletions", "Mutations", "Fusions", NA)) {
    event.type <- match.arg(event.type)
    if (is.null(events.mat)) {
        print(paste("Null", event.type, "matrix, skipping.."))
        return(NULL)
    }
    if (dim(events.mat)[1] == 0 | dim(events.mat)[2] == 0) {
        print(paste0("Not enough ", event.type,", skipping"))
        return(NULL)
    }
    
    # subset the two matrices to common samples
    common.samples <- intersect(colnames(vipermat), colnames(events.mat))
    vipermat <- vipermat[, common.samples]
    events.mat <- events.mat[, common.samples]
    
    # include only whitelist items
    if (!(is.na(whitelist))) {
        events.mat <- events.mat[intersect(rownames(events.mat), whitelist), ]
    }
    # filter to minmum number of somatic events
    if (!is.na(min.events)) {
        events.mat <- events.mat[apply(events.mat, 1, sum, na.rm = TRUE) >= min.events, ]
    }
    # test again after removing low freuquency events
    if (is.null(dim(events.mat))) {
        print(paste0("Not enough ", event.type,", skipping"))
        return(NULL)
    } else if (dim(events.mat)[1] == 0 | dim(events.mat)[2] == 0) {
        print(paste0("Not enough ", event.type,", skipping"))
        return(NULL)
    }
    
    nes <- areaEnrich(events.mat, vipermat, event.type)
    nes
}


#' aREA.enrich Compute aREA enrichment between all pairwise combinations of VIPER proteins and gene-level events
#' 
#' @param events.mat A Binary 0/1 matrix with columns as samples, and rows as proteins
#' @param vipermat A VIPER network of inferred activity scores with columns as samples, and rows as proteins
#' @param event.type Name of the event type for printing purposes
#' @return A matrix of enrichment scores with rows as event/gene names and columns as VIPER protein names
areaEnrich <- function(events.mat, vipermat, event.type) {
    
    # Convert mutations into a regulon-like object
    events.regulon <- apply(events.mat, 1, function(x) {
        l <- names(x[which(x == TRUE)])
        return(l)
    })
    
    # Calculate raw enrichment scores: each mutation against each TF columns are TFs, rownames are mutations
    es <- rea(t(vipermat), events.regulon, event.type = event.type)
    # Analytical correction
    dnull <- reaNULL(events.regulon)
    
    # Calculate pvalue of ES
    pval <- t(vapply(seq_len(length(dnull)), function(i, es, dnull) {
        dnull[[i]](es[i, ])$p.value
    }, es = es$groups, dnull = dnull, FUN.VALUE = numeric(ncol(es$groups))))
    
    # Convert the pvalues into Normalized Enrichment Scores
    nes <- qnorm(pval/2, lower.tail = FALSE) * es$ss
    # print (dim(nes)) print (length(events.regulon))
    rownames(nes) <- names(events.regulon)
    colnames(nes) <- rownames(vipermat)
    dimnames(pval) <- dimnames(nes)
    nes[is.na(nes)] <- 0
    # columns are TFs, rows are genomic events
    nes
}


#' This function calculates an Enrichment Score of Association based on how the features rank on the samples sorted by a specific gene
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @param eset Numerical matrix
#' @param regulon A list with genomic features as its names and samples as its entries, indicating presence of event
#' @param minsize The minimum number of events to use when calculating enrichment
#' @param maxsize The maximum number of events to use when calculating enrichment 
#' @param event.type Type of event being analyzed
#' @return A list containing two elements:
#' \describe{
#' \item{groups}{Regulon-specific NULL model containing the enrichment scores}
#' \item{ss}{Direction of the regulon-specific NULL model}
#' }
rea <- function(eset, regulon, minsize = 1, maxsize = Inf, event.type = NA) {
    # Filter for minimum sizes
    sizes <- vapply(regulon, length, FUN.VALUE = numeric(1))
    regulon <- regulon[sizes >= minsize]
    sizes <- vapply(regulon, length, FUN.VALUE = numeric(1))
    regulon <- regulon[sizes <= maxsize]
    
    temp <- unique(unlist(regulon))
    tw <- rep(1, length(temp))
    names(tw) <- temp
    
    # Calculations
    t <- eset
    t2 <- apply(t, 2, rank)/(nrow(t) + 1) * 2 - 1
    t1 <- abs(t2) * 2 - 1
    t1[t1 == (-1)] <- 1 - (1/length(t1))
    
    
    # tnorm
    t1 <- qnorm(t1/2 + 0.5)
    t2 <- qnorm(t2/2 + 0.5)
    
    if(is.na(event.type)) {
        event.type = "genomic events"
    }
    
    # Print some progress bar
    message("\nComputing associations of ", length(regulon)," ", event.type," with ", ncol(eset), " regulators")
    message("Process started at ", date())
    pb <- txtProgressBar(max = length(regulon), style = 3)
    
    temp <- lapply(seq_len(length(regulon)), function(i, regulon, t1, t2, tw, pb) {
        hitsamples <- regulon[[i]]
        hitsamples <- intersect(hitsamples, rownames(t1))
        
        # Mumbo-jumbo
        pos <- match(hitsamples, rownames(t1))
        
        heretw <- tw[match(hitsamples, names(tw))]
        
        sum1 <- matrix(heretw, 1, length(hitsamples)) %*% t2[pos, ]
        ss <- sign(sum1)
        ss[ss == 0] <- 1
        setTxtProgressBar(pb, i)
        sum2 <- matrix(0 * heretw, 1, length(hitsamples)) %*% t1[pos, ]
        return(list(es = as.vector(abs(sum1) + sum2 * (sum2 > 0))/sum(heretw), ss = ss))
    }, regulon = regulon, pb = pb, t1 = t1, t2 = t2, tw = tw)
    names(temp) <- names(regulon)
    message("\nProcess ended at ", date())
    es <- t(vapply(temp, function(x) x$es, FUN.VALUE = numeric(length(temp[[1]][[1]]))))
    ss <- t(vapply(temp, function(x) x$ss, FUN.VALUE = numeric(length(temp[[1]][[1]]))))
    colnames(es) <- colnames(ss) <- colnames(eset)
    return(list(groups = es, ss = ss))
}


#' This function generates the NULL model function, which computes the normalized enrichment score and associated p-value
#' 
#' @param regulon A list with genomic features as its names and samples as its entries
#' @param minsize Minimum number of event (or size of regulon) to calculate the model with 
#' @param maxsize Maximum number of event (or size of regulon) to calculate the model with 
#' @return A list of functions to compute NES and p-value
reaNULL <- function(regulon, minsize = 1, maxsize = Inf) {
    # Filter for minimum sizes
    sizes <- vapply(regulon, length, FUN.VALUE = numeric(1))
    regulon <- regulon[sizes >= minsize]
    sizes <- vapply(regulon, length, FUN.VALUE = numeric(1))
    regulon <- regulon[sizes <= maxsize]
    # complete list of all genes in any regulon
    temp <- unique(unlist(regulon))
    # list of all genes, weighted by 1
    tw <- rep(1, length(temp))
    names(tw) <- temp
    lapply(regulon, function(x, tw) {
        ww <- tw[match(x, names(tw))]
        ww <- ww/max(ww)
        # now it's a constant?
        ww <- sqrt(sum(ww^2))
        return(function(x, alternative = "two.sided") {
            x <- x * ww
            p <- switch(pmatch(alternative, c("two.sided", "less", "greater")), pnorm(abs(x), lower.tail = FALSE) * 2, pnorm(x, lower.tail = TRUE), pnorm(x, 
                                                                                                                                                          lower.tail = FALSE))
            list(nes = x, p.value = p)
        })
    }, tw = tw)
}

