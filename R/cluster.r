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
#' @keywords internal
clusterRange <- function(dis, range = c(2, 100), step = 1, cores = 1, method = c("pam", "kmeans"), data = NULL) {
  
  
  #idx <- 1
  result <- parallel::mclapply(range[1]:range[2], function(k, dis) {
    
    solution <- NULL
    clustering <- NULL
    centers <- NULL
    silinfo <- NULL
    if (method == "pam") {
      solution <- cluster::pam(dis, k, diss = TRUE, cluster.only = FALSE)
      clustering <- solution$clustering
      centers <- solution$medoids
      silinfo <- solution$silinfo
    } else if (method == "kmeans") {
      data[which(is.na(data))] <- 0
      solution <- stats::kmeans(data, k)
      clustering <- solution$cluster
      centers <- solution$centers
    } else {
      stop("Unrecognized clustering method")
    }
    
    cr <- clusterReliability(clustering, 1/as.matrix(dis), method = "global")
    element.cr <- clusterReliability(clustering, 1/as.matrix(dis), method = "element")
    cluster.cr <- clusterReliability(clustering, 1/as.matrix(dis), method = "cluster")
    
    
    ret <- NULL
    list(centers = centers, clustering = clustering, k = k, silinfo = silinfo,
         reliability = as.numeric(cr), element.reliability = element.cr, 
         cluster.reliability = cluster.cr)
    
  }, dis = dis, mc.cores = cores)
  
  result <- setNames(result, paste0(seq(range[1], range[2]), "clusters"))
  
  # create fields that have all the cluster reliability scores and all the avg silhouette scores
  rel.scores <- vapply(result, function(x){ x[["reliability"]] }, FUN.VALUE = double(1))
  sil.scores <- vapply(result, function(x){ x[["silinfo"]][["avg.width"]] }, FUN.VALUE = double(1))
  
  result[["all.cluster.reliability"]] <- rel.scores
  result[["all.sil.avgs"]] <- sil.scores
  
  result
}

#' Cluster membership reliability estimated by enrichment analysis
#'
#' This function estimates the cluster membership reliability using aREA
#'
#' @param cluster Vector of cluster memberships or list of cluster memberships
#' @param similarity Similarity matrix
#' @param xlim Optional vector of 2 components indicating the limits for computing AUC
#' @param method Character string indicating the mthod to compute reliability, 
#' either by element, by cluster or global
#' @return Reliability score for each element
#' @keywords internal
clusterReliability <- function(cluster, similarity, xlim = NULL, 
                               method = c("element", "cluster", "global")) {
  method <- match.arg(method)
  if (!is.list(cluster)) 
    cluster <- list(cluster = cluster)
  switch(method, element = {
    res <- lapply(cluster, function(cluster, similarity) {
      reg <- split(names(cluster), cluster)
      tmp <- sREA(similarity, reg)
      tmp <- lapply(seq_len(nrow(tmp)), function(i, tmp, reg) {
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
    res <- lapply(seq_len(length(cluster)), function(i, cluster, rel, xlim) {
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

#' Simple one-tail rank based enrichment analysis sREA
#' (for cluster analysis)
#' 
#' This function performs simple 1-tail rank based enrichment analysis
#' 
#' @param signatures Numeric matrix of signatures
#' @param groups List containing the groups as vectors of sample names
#' @return Matrix of Normalized Enrichment Zcores
#' @keywords internal
sREA <- function(signatures, groups) {
  if (is.null(nrow(signatures))) 
    signatures <- matrix(signatures, length(signatures), 1, dimnames = list(names(signatures), "sample1"))
  # ranked signatures matrix: samples are rows. Rank
  sig <- qnorm(apply(signatures, 2, rank)/(nrow(signatures) + 1))
  gr <- vapply(groups, function(x, samp) {
    samp %in% x
  }, samp = rownames(sig), FUN.VALUE = logical(length(unlist(groups))))
  gr <- t(gr)
  # non negative counts
  nn <- rowSums(gr)
  # fractions of prev computation rows are groups
  gr <- gr/nn
  es <- gr %*% sig
  return(es * sqrt(nn))
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
#' @keywords internal
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
#' @keywords internal
integrateTZ <- function(x, y) {
  pos <- order(x)
  x <- x[pos]
  y <- y[pos]
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
}

#' Get best clustering supervised
#' 
#' This function calculates the best clustering based on optimizing for survival
#' separation in the event there are statistically equivalent solutions
#' 
#' @param search.results results object from clusterRange with all the information
#' for each possible clustering solution
#' @param top.sol the top solution from the first pass of the analysis
#' @param survival.data dataframe with survival data for each sample
#' @param progression.free T/F about whether to use progression free instead of overall survival
#' @return best clustering solution k
#' @keywords internal
getBestClusteringSupervised <- function(search.results, top.sol, survival.data, progression.free) {
  
  best.clustering <- search.results[[top.sol]]
  
  # Do KS test to compare distributions of sample reliability scores
  pvals <- vapply(2:10, function(k) {
    cs <- search.results[[(k-1)]]
    pval.1 <- suppressWarnings(ks.test(cs$element.reliability, best.clustering$element.reliability))$p.value
    pval.1
  }, FUN.VALUE = numeric(1))
  
  names(pvals) <- 2:10
  fwers <- sort(p.adjust(pvals, method='bonferroni'))
  # those where we can't reject the null hypothesis are considered equivalent
  equivalent.clusters <- as.numeric(names(fwers[which(fwers > 0.05)]))

  # TODO: presumably if there is only one equivalent cluster, it's the original one
  # can skip if so
  if(length(equivalent.clusters) == 1) {
    message("No statistically equivalent clustering solutions. Keeping this k:", equivalent.clusters)
    return(equivalent.clusters)
  }
  
  to.print <- paste0(equivalent.clusters, collapse = " , ")
  message("Clustering solution(s) [" , to.print, "] are statistically equivalent to the top solution. Now checking for best survival")
  
  
  # Confirm survival data has matching sample names 
  name.res <- checkSurvivalNames(search.results, survival.data)
  
  # For each equivalent clustering solution do survival analysis to find best and worst
  
  new.best.clustering <- bestVSworst(equivalent.clusters, search.results, survival.data, name.res, progression.free)
  
  new.best.clustering
  
}

#' Check survival names
#' 
#' Function to check names of samples in survival data 
#' Needed mainly because names are shorter in patient data
#' @param search.results clustering object from clusterRange
#' @param survival.data dataframe with survival data for each sample
#' @return list with length of sample names and column name for samples
#' @keywords internal
checkSurvivalNames <- function(search.results, survival.data) {
  
  clus.sample.names <- names(search.results[[1]]$clustering)
  
  cname <- intersect(c("sample", "patient", "bcr_patient_barcode"), colnames(survival.data))
  if(length(cname) == 0) {
    stop("Survival data does not have correct column name for samples. Rename to either 'sample' or 'patient'")
  }
  
  survival.sample.names <- unlist(survival.data[,cname], use.names = F)
  overlap <- intersect(clus.sample.names, survival.sample.names)
  
  # if there is some overlap the names must be in the same format so continue
  if(length(overlap) > 0) {
    return(list(len = nchar(overlap[1]), cname = cname))
  } 
  
  # if no overlap check if one set of names is a subset of the other
  smallest.len <- min(nchar(clus.sample.names[1]), nchar(survival.sample.names[1]))
  
  clus.sample.names <- sampleNameFilter(clus.sample.names, smallest.len)
  survival.sample.names <- sampleNameFilter(survival.sample.names, smallest.len)
  
  overlap <- intersect(clus.sample.names, survival.sample.names)
  
  # if there is some overlap save the length that the names need to be in order to match
  if(length(overlap) > 0) {
    return(list(len = nchar(overlap[1]), cname = cname))
  } else {
    stop("Sample names do not match. Cannot do survival analysis")
  }
  
  length(list(len = nchar(overlap[1]), cname = cname))
  
}


#' Best vs Worst 
#' 
#' Take clustering solutions with equivalent reliability distributions to the best
#' one and calculate which result has the highest p value between the best and 
#' worst surviving clusters
#' @param equivalent.clusters vector indicating k of equivalent clusters
#' @param search.results clustering object from clusterRange
#' @param survival.data dataframe with survival data for each sample
#' @param name.res list with length of sample names and column name for samples
#' @param progression.free T/F about whether to use progression free instead of overall survival
#' @return k for the top solution
#' @keywords internal
#' @import survival
bestVSworst <- function(equivalent.clusters, search.results, survival.data, name.res, progression.free) {
  
  
  res.pvals <- vapply(equivalent.clusters, function(k) {
    # get clustering solution and join to the survival dataframe
    # confirm number of samples that have survival data
    clustering <- search.results[[(k-1)]]$clustering
    
    names(clustering) <- sampleNameFilter(names(clustering), name.res$len)
    clustering <- tibble::enframe(clustering, "samples", "cluster")
    
    num.samples <- nrow(clustering)
    
    join.cols <- setNames(name.res$cname, "samples")
    full.df <- dplyr::inner_join(clustering, survival.data, by = join.cols)
    
    num.w.surv <- nrow(full.df)
    
    message("Survival data found for ", num.w.surv, " out of ", num.samples, " samples.")
    
    if(num.w.surv == 0) {
      stop("Survival data not properly added to clusters. Check inputs.")
    }
    
    # Make survival object
    # Check for "time" and "event" as column names in the survival data
    # if not look for OS/OS.time or PFI/PFI.time
    
    full.df <- survfit(full.df, progression.free)
    
    # Find the best and worst surviving clusters based on the statistics
    SurvDiff <- survival::survdiff(survObj ~ cluster, data=full.df, rho=0)
    survival.scores <- SurvDiff$obs/SurvDiff$exp
    best.surv.cluster <- which.min(survival.scores)
    worst.surv.cluster <- which.max(survival.scores)
    
    # Redo the survival analysis difference between just these clusters
    clin.subset <- full.df %>% dplyr::filter(cluster %in% c(best.surv.cluster, worst.surv.cluster))
    clin.subset <- survfit(clin.subset, progression.free)
    SurvDiff <- survival::survdiff(survObj ~ cluster, data=clin.subset, rho=0)
    pval <- 1 - pchisq(SurvDiff$chisq, length(SurvDiff$n) - 1)
    pval
    
  }, FUN.VALUE = numeric(1)) 
  
  names(res.pvals) <- equivalent.clusters
  
  best.k <- names(which.min(res.pvals))
  
  message("Optimal clustering solution based on survival separation is: ", best.k)
  
  best.k
  
}

#' Survfit 
#' 
#' Helper function to attach survObj to clinical data 
#' @param full.df dataframe of survival data with clusters added
#' @param progression.free whether or not to do overall survival or progression free
#' @return dataframe with survObj attached
#' @import survival
#' @keywords internal
survfit <- function(full.df, progression.free) {
  
  if(all(c("time", "event") %in% colnames(full.df))) {
    full.df$survObj <- survival::Surv(full.df$time, full.df$event)
  } else if (isFALSE(progression.free) & all(c("OS.time", "OS") %in% colnames(full.df))) {
    full.df$survObj <- survival::Surv(full.df$OS.time, full.df$OS)
  } else if (isTRUE(progression.free) & all(c("PFI.time", "PFI") %in% colnames(full.df))) {
    full.df$survObj <- survival::Surv(full.df$PFI.time, full.df$PFI)
  } else {
    stop("Survival data not properly formatted. Make sure the column names are 'time' and 'event'.")
  }
  
  full.df
  
}







