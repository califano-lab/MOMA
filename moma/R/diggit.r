
library(qvalue)
library(cluster)

#' @title Use Stouffer's method to combine z-scores of DIGGIT interactions for each cMR protein. Combines only positively associated DIGGIT scores by default.  
#' @import stats
#' @param interactions A list indexed by TF, includes z-scores or p-values for each interacting event
#' @param from.p Integrate p-values or z-scores (default z-scores; from.p = FALSE)
#' @param pos.nes.only Use only positive NES scores to rank proteins (default TRUE)
#' @return A list indexed by TF, a stouffer integrated z-score
stouffer.integrate.diggit <- function(interactions, from.p = FALSE, pos.nes.only = TRUE) {
    
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


#' @title Filter interactions from NES (DIGGIT) scores and corresponding background-corrected scores. Use this version in the Bayes model to rank TFs
#' @import stats
#' @param corrected.scores A list indexed by the genomic event/gene with corresponding pvals and qvals for
#' each TF
#' @param nes.scores Matrix with tfs as columns, rows are genomic events
#' @param cindy CINDy algorithm output matrix
#' @param p.thresh P-value threshold (default=0.05)
#' @param cindy.only Consider only CINDy validated interactions (default=TRUE)
#' @return a list (indexed by VIPER protein) of significant genomic interactions 
#' and associated pvals over the background (null TF) model, and NES scores
sig.interactors.DIGGIT <- function(corrected.scores, nes.scores, cindy, p.thresh = 0.05, cindy.only = TRUE) {
    
    pvals.matrix <- get.pvals.matrix(corrected.scores)
    
    # input validation
    if (!is.numeric(p.thresh)) {
        print("Error: invalid value supplied for p-value threshold!")
        q()
    }
    
    ## Apply joint NES + p-value over background (null TF) threshold over each Viper Protein return the raw NES scores only for those significant over the
    ## background and including CINDy (if applicable)
    viper.interactors <- lapply(colnames(pvals.matrix), function(viperProt) {
        
        # find the over-null-TF-background scores with an significant, uncorrected p-value print (paste('Processing : ', viperProt))
        pvals <- as.numeric(pvals.matrix[, as.character(viperProt)])
        nes.vec <- as.numeric(nes.scores[, as.character(viperProt)])
        # print (nes.vec) subset to significant p-values
        row.idx <- which(pvals < p.thresh)
        pvals <- pvals[row.idx]
        names(pvals) <- rownames(pvals.matrix)[row.idx]
        
        # if (!all(names(pvals) == names(nes.vec))) { print ('Error: data not aligned for aREA / aREA corrected p-values') }
        
        # subset the NES vector for this TF, threshold again on NES scores as a sanity check on the basic enrichment (i.e. remove those with high
        # over-background scores simply because the background is de-enriched)
        nes.vec <- nes.scores[, which(colnames(nes.scores) == as.character(viperProt))]
        nes.vec <- nes.vec[which(2 * (1 - pnorm(abs(nes.vec))) < p.thresh)]
        
        fusion.index <- unlist(lapply(names(nes.vec), function(x) if (length(strsplit(x, "_")[[1]]) > 1) 
            TRUE else FALSE))
        # subset to CINDY validated upstream regulators
        if (cindy.only && is.null(cindy[[viperProt]])) {
            return(c())
        } else if (cindy.only) {
            
            # lookup takes about a second...
            cindy.thisTF <- cindy[[viperProt]]
            
            upstream.cindy.modulators <- names(cindy.thisTF)
            cindy.nes.vec <- na.omit(nes.vec[match(upstream.cindy.modulators, names(nes.vec))])
            if (length(cindy.nes.vec) == 0) {
                return(c())
            }
            
            # add fusions in, they don't have CINDy scores
            fus.vec <- nes.vec[fusion.index]
            entrez.vec <- cindy.nes.vec[which(!is.na(cindy.nes.vec))]
            
            # Independence of these: multiply through
            cindy.entrez.pvals <- cindy.thisTF[(names(cindy.thisTF) %in% names(entrez.vec))]
            entrez.pvals <- 2 * (1 - pnorm(abs(entrez.vec)))
            entrez.pvals <- entrez.pvals * cindy.entrez.pvals
            
            entrez.vec.corrected <- (1 - qnorm(entrez.pvals)) * sign(entrez.vec)
            nes.vec <- c(entrez.vec.corrected, fus.vec)
        }
        
        # keep interactions significant above the null model and with a significant raw aREA score
        nes.vec <- nes.vec[intersect(names(nes.vec), names(pvals))]
        nes.vec <- nes.vec[which(!is.na(nes.vec))]
        
        # print (paste('num interactions:', length(nes.vec)))
        nes.vec <- sort(nes.vec, decreasing = TRUE)
        nes.vec
    })
    names(viper.interactors) <- colnames(pvals.matrix)
    viper.interactors
}


#' @title Use 'aREA' to calculate the enrichment between each genomic event - VIPER inferred protein pair. Requires pre-computed VIPER scores and a binary events matrix. Will use only samples in both event and VIPER matrices. 
#' @param vipermat Pre-computed VIPER scores with samples as columns and proteins as rows 
#' @param events.mat Binary 0/1 events matrix with samples as columns and genes or events as rows
#' @param whitelist Only compute associations for events in this list
#' @param blacklist Exclude associations for events in this list
#' @param min.events Only compute enrichment if the number of samples with these events is GTE to this
#' @return A matrix of aREA scores, dimensions are nrow(events.mat) x nrow(vipermat) 
associate.events <- function(vipermat, events.mat, min.events = NA, whitelist = NA, blacklist = NA) {
    
    if (is.null(events.mat)) {
        print("Null mutation matrix, skipping..")
        return(NULL)
    }
    if (dim(events.mat)[1] == 0 | dim(events.mat)[2] == 0) {
        print("Not enough mutations...skipping")
        return(NULL)
    }
    
    # subset the two matrices to common samples
    common.samples <- intersect(colnames(vipermat), colnames(events.mat))
    vipermat <- vipermat[, common.samples]
    events.mat <- events.mat[, common.samples]
    
    # remove blacklist
    if (all(is.na(blacklist))) {
        # print (blacklist)
        events.mat <- events.mat[setdiff(rownames(events.mat), blacklist), ]
    }
    # include only whitelist items
    if (!all(is.na(whitelist))) {
        events.mat <- events.mat[intersect(rownames(events.mat), whitelist), ]
    }
    # filter to minmum number of somatic events
    if (!is.na(min.events)) {
        events.mat <- events.mat[apply(events.mat, 1, sum, na.rm = TRUE) >= min.events, ]
    }
    # test again after removing low freuquency events
    if (is.null(dim(events.mat))) {
        print("Not enough events...skipping")
        return(NULL)
    } else if (dim(events.mat)[1] == 0 | dim(events.mat)[2] == 0) {
        print("Not enough events...skipping")
        return(NULL)
    }
    
    nes <- aREA.enrich(events.mat, vipermat)
    nes
}

#' @title Compute the empirical q-values of each genomic-event/VIPER gene pair against the background distribution of associations with a given set of 'null' VIPER genes (i.e. low activity TFs)
#' @param vipermat viper inferences matrix, samples are columns, rows are TF entrez gene IDs 
#' @param nes scores for each mutation (rows) against each TF (columns) 
#' @param null.TFs low-importance TFs used to calculate null distributions
#' @param alternative Alternative defaults to 'both' : significant p-values can come from both sides of the null distribution 
#' @return A named list of qvalues for each TF/cMR protein. Each entry contains a vector of q-values for all associated events; names are gene ids 
get.diggit.empiricalQvalues <- function(vipermat, nes, null.TFs, alternative = "both") {
    
    # subset NES to Viper Proteins in the vipermat only
    nes <- nes[, as.character(rownames(vipermat))]
    
    nes.em.qvals <- apply(nes, 1, function(x, alternative) {
        
        null.VEC <- x[as.character(null.TFs)]
        null.VEC <- null.VEC[which(!is.na(null.VEC))]
        # get empirical q-values for both upper and lower tails of NES / DIGGIT statistics
        qvals <- get.empirical.qvals(x, null.VEC, alternative)
        qvals
    }, alternative = alternative)
    
    names(nes.em.qvals) <- rownames(nes)
    nes.em.qvals
}



#' @title Compute aREA enrichment for the proteins in a given regulon, against vipermat scores supplied
#' @param regulon ARACNE regulon object
#' @param vipermat A VIPER network of inferred activity scores with columns as patient samples, and rows as proteins
#' @return A matrix of enrichment scores with rows as event/gene names and columns as VIPER protein names
aREA.regulon_enrich <- function(regulon, vipermat) {
    
    regulon <- lapply(regulon, function(x) as.character(x))
    # Calculate raw enrichment scores: each mutation against each TF columns are TFs, rownames are mutations
    es <- rea(t(vipermat), regulon)
    # Analytical correction
    dnull <- reaNULL(regulon)
    
    # Calculate pvalue of ES
    pval <- t(vapply(1:length(dnull), function(i, es, dnull) {
        dnull[[i]](es[i, ])$p.value
    }, es = es$groups, dnull = dnull, FUN.VALUE = numeric(ncol(es$groups))))
    
    # Convert the pvalues into Normalized Enrichment Scores
    nes <- qnorm(pval/2, lower.tail = FALSE) * es$ss
    # print (dim(nes)) print (length(regulon))
    rownames(nes) <- names(regulon)
    colnames(nes) <- rownames(vipermat)
    dimnames(pval) <- dimnames(nes)
    nes[is.na(nes)] <- 0
    # columns are TFs, rows are genomic events
    nes
}


#' @title aREA.enrich Compute aREA enrichment between all pairwise combinations of VIPER proteins and gene-level events
#' @param events.mat A Binary 0/1 matrix with columns as samples, and rows as proteins
#' @param vipermat A VIPER network of inferred activity scores with columns as samples, and rows as proteins
#' @return A matrix of enrichment scores with rows as event/gene names and columns as VIPER protein names
#' @examples 
#' nes <- aREA.enrich(moma$gbm.example[['rawsnp']], moma$gbm.example[['vipermat']])
#' dim(nes)  
#' @export
aREA.enrich <- function(events.mat, vipermat) {
    
    # Convert mutations into a regulon-like object
    events.regulon <- apply(events.mat, 1, function(x) {
        l <- names(x[which(x == TRUE)])
        return(l)
    })
    
    # Calculate raw enrichment scores: each mutation against each TF columns are TFs, rownames are mutations
    es <- rea(t(vipermat), events.regulon)
    # Analytical correction
    dnull <- reaNULL(events.regulon)
    
    # Calculate pvalue of ES
    pval <- t(vapply(1:length(dnull), function(i, es, dnull) {
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


#' @title This function calculates an Enrichment Score of Association based on how the features rank on the samples sorted by a specific gene
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @param eset Numerical matrix
#' @param regulon A list with genomic features as its names and samples as its entries, indicating presence of event
#' @param minsize The minimum number of events to use when calculating enrichment
#' @param maxsize The maximum number of events to use when calculating enrichment 
#' @return A list containing two elements:
#' \describe{
#' \item{groups}{Regulon-specific NULL model containing the enrichment scores}
#' \item{ss}{Direction of the regulon-specific NULL model}
#' }
rea <- function(eset, regulon, minsize = 1, maxsize = Inf) {
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
    
    # Print some progress bar
    message("\nComputing associations of ", length(regulon), " genomic events with ", ncol(eset), " genes")
    message("Process started at ", date())
    pb <- txtProgressBar(max = length(regulon), style = 3)
    
    temp <- lapply(1:length(regulon), function(i, regulon, t1, t2, tw, pb) {
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



#' @title This function generates the NULL model function, which computes the normalized enrichment score and associated p-value
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


#' @title Utility function
#' @param corrected.scores - corrected p-values processed by 'qvals' package
#' @return A matrix of p-values for scores between genes/events (rows) and TFs (columns)
get.pvals.matrix <- function(corrected.scores) {
    # order of VIPER proteins/TFs
    tf.names.order <- names(corrected.scores[[1]]$qvals)
    pvals.matrix <- matrix(unlist(lapply(corrected.scores, function(x) {
        pvals <- x$pvals[tf.names.order]
        pvals
    })), byrow = TRUE, ncol = length(tf.names.order))
    
    colnames(pvals.matrix) <- tf.names.order
    rownames(pvals.matrix) <- names(corrected.scores)
    pvals.matrix
}

#' @title Utility function
#' @param vipermat - matrix of VIPER scores with columns as samples, rows as protein names
#' @param fdr.thresh - BH-FDR threshold (default 0.05 FDR rate)
#' @return A vector of normalized z-scores, named by TF id
viper.getTFScores <- function(vipermat, fdr.thresh = 0.05) {
    
    # for each gene, count the number samples with scores for each, and weight by that
    w.counts <- apply(vipermat, 1, function(x) {
        data.counts <- length(which(!is.na(x)))
        data.counts
    })
    w.counts <- w.counts/ncol(vipermat)
    
    vipermat[is.na(vipermat)] <- 0
    
    # normalize element scores to sum to 1 (optional - use weighted element scores based on silhouette)
    element.scores <- rep(1, ncol(vipermat))
    element.scores <- element.scores/sum(element.scores)
    
    # mean weighted VIPER score across samples
    w.means <- apply(vipermat, 1, function(x) {
        res <- sum(x * element.scores)
        res
    })
    # weight by the counts for each
    w.means <- w.means * w.counts
    names(w.means) <- rownames(vipermat)
    
    # only look at those with positive (high) score w.means <- sort(w.means[which(w.means > 0)], decreasing=TRUE)
    
    zscores <- w.means
    zscores
}

#' @title Calculate p-values from pseudo zscores / VIPER aREA scores, threshold
#' @param zscores Vector of normally distributed z-scores representing protein activities. 
#' @param fdr.thresh Threshold for false discovery rate, default is 0.05
#' @return Get the names of proteins with significant z-scores, after multi-hypothesis correction
viper.getSigTFS <- function(zscores, fdr.thresh = 0.05) {
    
    # calculate pseudo-pvalues and look at just significant pvals/scores
    pvals <- -pnorm(abs(zscores), log.p = TRUE) * 2
    pvals[which(pvals > 1)] <- 1
    # correct unless option is NULL
    sig.idx <- which(p.adjust(pvals, method = "BH") < fdr.thresh)
    pvals <- pvals[sig.idx]
    
    names(pvals)
}

#' @title Retain TCGA sample ids without the final letter designation ('A/B/C')
#' @param mat Matrix of expression or protein activity scores. Columns are sample names, rows are genes
#' @return An identical matrix with new (shorter) column names
samplename.filter <- function(mat) {
    # filter down to sample Id without the 'A/B/C sample class'.
    sample.ids <- vapply(colnames(mat), function(x) substr(x, 1, 15), FUN.VALUE = character(1))
    colnames(mat) <- sample.ids
    mat
}

#' @title get.empirical.qvals
#' @param test.statistics P-values generated from the test comparisons
#' @param null.statistics P-values generated under the null (permutation) model
#' @param alternative Optional : 1 or 2 tails used to generate the p-value (default='both')
#' @return A list with both the qvalues and empirical p-values from the supplied test and null stats
get.empirical.qvals <- function(test.statistics, null.statistics, alternative = "both") {
    
    # calculate the upper and lower tail
    if (alternative == "both") {
        
        test.statistics <- sort(abs(test.statistics), decreasing = TRUE)
        null.statistics <- abs(null.statistics)
        
        em.pvals <- qvalue::empPvals(test.statistics, null.statistics)
        qvals <- rep(1, length(em.pvals))
        tryCatch({
            qvals <- qvalue::qvalue(em.pvals)$qvalue
        }, error = function(e) {
            # if pi0, the estimated proportion of true null hypothesis <= 0, it might fail: in that case set to zero and return p-values anyways
            qvals <- rep(1, length(em.pvals))
        })
        names(qvals) <- names(test.statistics)
        names(em.pvals) <- names(test.statistics)
        return(list(qvals = qvals, pvals = em.pvals))
    } else {
        stop(paste(" alternative ", alternative, " not implemented yet!"))
    }
}

#' Simple one-tail rank based enrichment analysis sREA
#' 
#' This function performs simple 1-tail rank based enrichment analysis
#' 
#' @param signatures Numeric matrix of signatures
#' @param groups List containing the groups as vectors of sample names
#' @return Matrix of Normalized Enrichment Zcores
sREA <- function(signatures, groups) {
    if (is.null(nrow(signatures))) 
        signatures <- matrix(signatures, length(signatures), 1, dimnames = list(names(signatures), "sample1"))
    # ranked signatures matrix: samples are rows. Rank
    sig <- qnorm(apply(signatures, 2, rank)/(nrow(signatures) + 1))
    gr <- sapply(groups, function(x, samp) {
        samp %in% x
    }, samp = rownames(sig))
    gr <- t(gr)
    # non negative counts
    nn <- rowSums(gr)
    # fractions of prev computation rows are groups
    gr <- gr/nn
    es <- gr %*% sig
    return(es * sqrt(nn))
}

#' Cluster membership reliability estimated by enrichment analysis
#'
#'This function etimates the cluster membership reliability using aREA
#'
#' @param cluster Vector of cluster memberships or list of cluster memberships
#' @param similarity Similarity matrix
#' @param xlim Optional vector of 2 components indicating the limits for computing AUC
#' @param method Character string indicating the mthod to compute reliability, either by element, by cluster or global
#' @return Reliability score for each element
clusterReliability <- function(cluster, similarity, xlim = NULL, method = c("element", "cluster", "global")) {
    method <- match.arg(method)
    if (!is.list(cluster)) 
        cluster <- list(cluster = cluster)
    switch(method, element = {
        res <- lapply(cluster, function(cluster, similarity) {
            reg <- split(names(cluster), cluster)
            tmp <- sREA(similarity, reg)
            tmp <- lapply(1:nrow(tmp), function(i, tmp, reg) {
                tmp[i, ][colnames(tmp) %in% reg[[which(names(reg) == (rownames(tmp)[i]))]]]
            }, tmp = tmp, reg = reg)
            res <- unlist(tmp, use.names = FALSE)
            names(res) <- unlist(lapply(tmp, names), use.names = FALSE)
            return(res[match(names(cluster), names(res))])
        }, similarity = similarity)
        if (length(res) == 1) return(res[[1]])
        return(res)
    }, cluster = {
        rel <- clusterReliability(cluster, similarity, method = "element")
        if (!is.list(rel)) rel <- list(cluster = rel)
        if (is.null(xlim)) xlim <- range(unlist(rel, use.names = FALSE))
        res <- lapply(1:length(cluster), function(i, cluster, rel, xlim) {
            tapply(rel[[i]], cluster[[i]], function(x, xlim) {
                1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps = 1000)/diff(xlim)
            }, xlim = xlim)
        }, cluster = cluster, rel = rel, xlim = xlim)
        if (length(res) == 1) return(res[[1]])
        return(res)
    }, global = {
        rel <- clusterReliability(cluster, similarity, method = "element")
        if (!is.list(rel)) rel <- list(cluster = rel)
        if (is.null(xlim)) xlim <- range(unlist(rel, use.names = FALSE))
        res <- lapply(rel, function(x, xlim) {
            1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps = 1000)/diff(xlim)
        }, xlim = xlim)
        if (length(res) == 1) return(res[[1]])
        return(res)
    })
}

#' Cluster Range
#' 
#' This function generate an cluster structure with 'k' groups and computes the cluster reliability score
#' where 'k' is a range of values 
#' 
#' @param dis Distance object
#' @param range vector with start and end 'k' 
#' @param step Integer indicating the incremental number of clusters to add in each iteration
#' @param cores Maximum number of CPU cores to use
#' @param method Either 'pam' k-mediods or kmeans. Must supply the original data matrix if using kmeans
#' @param data Original data matrix
#' @return list of cluster reliability scores by 'k', 'clustering' (the vector solution) and 'reliability' 
#' as well as 'medoids' labels
clusterRange <- function(dis, range = c(1, 100), step = 1, cores = 1, method = c("pam", "kmeans"), data = NULL) {
    
    debug = TRUE
    idx <- 1
    result <- mclapply(range[1]:range[2], function(k, dis) {
        
        solution <- NULL
        clustering <- NULL
        centers <- NULL
        if (method == "pam") {
            solution <- pam(dis, k, diss = TRUE, cluster.only = FALSE)
            clustering <- solution$clustering
            centers <- solution$medoids
        } else if (method == "kmeans") {
            data[which(is.na(data))] <- 0
            solution <- kmeans(data, k)
            clustering <- solution$cluster
            centers <- solution$centers
        } else {
            print("Error: unrecognized clustering method")
            return(NULL)
        }
        
        cr <- clusterReliability(clustering, 1/as.matrix(dis), method = "global")
        element.cr <- clusterReliability(clustering, 1/as.matrix(dis), method = "element")
        cluster.cr <- clusterReliability(clustering, 1/as.matrix(dis), method = "cluster")
        if (debug) {
            print(paste0(k, " ", as.numeric(cr)))
        }
        
        ret <- NULL
        list(centers = centers, clustering = clustering, k = k, reliability = as.numeric(cr), element.reliability = element.cr, cluster.reliability = cluster.cr)
        
    }, dis = dis, mc.cores = cores)
    
    result
}

#' Numerical integration of functions
#' 
#' Integrates numerically a function over a range using the trapezoid method
#' 
#' @param f Function of 1 variable (first argument)
#' @param xmin Number indicating the min x value
#' @param xmax Number indicating the max x value
#' @param steps Integer indicating the number of steps to evaluate
#' @param ... Additional arguments for \code{f}
#' @return Number
integrateFunction <- function(f, xmin, xmax, steps = 100, ...) {
    x <- seq(xmin, xmax, length = steps)
    y <- f(x, ...)
    integrateTZ(x, y)
}

#' Integration with trapezoid method
#' 
#' This function integrate over a numerical range using the trapezoid method
#' 
#' @param x Numeric vector of x values
#' @param y Numeric vector of y values
#' @return Number
integrateTZ <- function(x, y) {
    pos <- order(x)
    x <- x[pos]
    y <- y[pos]
    idx = 2:length(x)
    return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
}


