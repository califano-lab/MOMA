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
clusterRange <- function(dis, range = c(2, 100), step = 1, cores = 1, method = c("pam", "kmeans"), data = NULL) {
  
  #set.seed(1, kind = "L'Ecuyer-CMRG")
  #debug = TRUE
  idx <- 1
  result <- parallel::mclapply(range[1]:range[2], function(k, dis) {
    
    solution <- NULL
    clustering <- NULL
    centers <- NULL
    if (method == "pam") {
      solution <- cluster::pam(dis, k, diss = TRUE, cluster.only = FALSE)
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
    #if (debug) {
    #   print(paste0(k, " ", as.numeric(cr)))
    #}
    
    ret <- NULL
    list(centers = centers, clustering = clustering, k = k, reliability = as.numeric(cr), element.reliability = element.cr, cluster.reliability = cluster.cr)
    
  }, dis = dis, mc.cores = cores)
  
  result
}

#' Cluster membership reliability estimated by enrichment analysis
#'
#' This function estimates the cluster membership reliability using aREA
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

