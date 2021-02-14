#' Sample Event Overlap
#' 
#' Updated function to compute which sample-specific alterations overlap with 
#' genomic events that are explained via DIGGIT/upstream of top regulators.
#' @param momaObj main object holding all the previously calculated results
#' @param viper.samples vector of samples to be analyzed
#' @param selected.tfs top ranking tfs that are active in the samples under consideration
#' @param interaction.df dataframe of all interactions associated with the regulators
#' @param cytoband.collapse T/F of whether to consolidate cnvs in the same cytoband to be considered as one event
#' @param topN how many regulators to consider for the saturation analysis 
#' @param verbose T/F whether to print sample names
#' @return a list obj with all the coverage information for each sample analyzed
#' @keywords internal
sampleEventOverlap <- function(momaObj, viper.samples, selected.tfs,
                               interaction.df, cytoband.collapse = T, topN = 100, verbose = F) {

  # restrict selected.tfs to topN, starting with 100 may change later
  selected.tfs <- selected.tfs[seq_len(topN)]
  
  # confirm DIGGIT has been run already
  if (is.null(momaObj$hypotheses)) {
    stop("Error: no hypothesis set for momaRunner class object!!")
  }
  
  event.types <- names(momaObj$hypotheses$events)
  
  # consider only samples triple intersection of samples with CNV, mutation and RNA (i.e. VIPER inferences) data
  samples.w.genomics <- intersect(colnames(momaObj$mut), colnames(momaObj$cnv))
  all.sample.names <- intersect(samples.w.genomics, viper.samples)
  
  
  
  coverage <- lapply(all.sample.names, function(sample) {
    
    if(isTRUE(verbose)){
      message("Computing coverage for sample ", sample)
    } 
    
    # find active and inactive proteins in this sample
    viper.pvals <- 2*pnorm(-abs(momaObj$viper[, sample]))
    active.proteins <- names(viper.pvals[viper.pvals < 0.05]) 
    if (length(active.proteins) == 0) {
      warning(paste0("No active proteins found for sample: ", sample))
    }
    
    # collect all genomic events for this sample
    # for amp/dels keep cytoband location for consolidation later
    
    all.sample.genomics <- tibble::tibble(.rows = 0)
    for (etype in event.types) {
      
      # First double check that sample exists in data set. Mainly for fusions.
      # Could consider making this the only check and allowing samples with  
      # partial genomic inputs to be considered...
      if(!sample %in% colnames(momaObj$hypotheses$matrices[[etype]])) next
      
      present <- as.character(names(which(momaObj$hypotheses$matrices[[etype]][,sample] > 0)))
      
      df <- tibble::tibble(event = present, type = etype)
      all.sample.genomics <- dplyr::bind_rows(all.sample.genomics, df)
    }
   
    # TODO: Fix in case of fusions....? 
    all.sample.genomics <- dplyr::left_join(all.sample.genomics, momaObj$gene.loc.mapping, by = c("event" = "Entrez.IDs")) 
    
    ## total number of events +/- cytoband collapse
    if(isTRUE(cytoband.collapse)) {
      cnv.events <- all.sample.genomics %>% 
        dplyr::filter(.data$type %in% c("amp", "del")) %>% 
        dplyr::select(.data$Cytoband) %>% unique() %>% nrow()
      non.cnv.events <- all.sample.genomics %>% 
        dplyr::filter(!.data$type %in% c("amp", "del")) %>% 
        nrow()
      total.events <- cnv.events + non.cnv.events
    } else {
      total.events <- nrow(all.sample.genomics)
    }
    
    ## Sample coverage calculation
    ## Going to start out looking through top 100 regulators, may change later
    
    sample.cover <- lapply(seq_along(selected.tfs), function(x) c())
    idx.range <- seq_along(selected.tfs)
    
    # store the set of active MRs, if it's identical to the last iteration, update and skip re-computation
    active.mrs.length.LAST <- 0
    last.idx <- 1
    
    for(i in idx.range) {
      # use these for coverage
      top.N.tfs <- selected.tfs[seq_len(i)]
      active.mrs <- intersect(active.proteins, top.N.tfs)
      
      # if no extra active mrs this round, copy previous stats and continue
      if ((length(active.mrs) == active.mrs.length.LAST) & (i > 1)) {
        sample.cover[[i]] <- sample.cover[[last.idx]]
        next
      } else {
        active.mrs.length.LAST <- length(active.mrs)
        last.idx <- i
      }
      
      # get events associated with this regulator from the interactions.df see which are in this sample
      active.mrs.interactions <- interaction.df %>% dplyr::filter(.data$regulator %in% active.mrs)
      
      covered.sample.genomics <- tibble::tibble(.rows = 0)
      for(etype in event.types) {
        sample.events <- all.sample.genomics %>% 
          dplyr::filter(.data$type == etype) %>%
          dplyr::select(.data$event) %>% 
          unlist(use.names = F) %>% unique()
        df <- active.mrs.interactions %>% 
          dplyr::filter(.data$type == etype & .data$event %in% sample.events) 
        
        covered.sample.genomics <- dplyr::bind_rows(covered.sample.genomics, df) %>% 
          dplyr::distinct()
      }
      
      # count events taken explained
      # consolidate cnvs HERE if necessary to match up with how original count was taken
      
      events.explained <- covered.sample.genomics %>% 
        dplyr::select(.data$event, .data$type, .data$Cytoband) %>% 
        dplyr::distinct()
      
      if(isTRUE(cytoband.collapse)) {
        cnv.events <- events.explained %>% 
          dplyr::filter(.data$type %in% c("amp", "del")) %>% 
          dplyr::select(.data$Cytoband, .data$type) %>% 
          dplyr::rename(event = .data$Cytoband) %>%
          unique() 
        non.cnv.events <- events.explained %>% 
          dplyr::filter(!.data$type %in% c("amp", "del")) %>% 
          dplyr::select(.data$event, .data$type) %>%
          unique()
        
        events.explained <- dplyr::bind_rows(cnv.events, non.cnv.events)
          
        event.count <- nrow(cnv.events) + nrow(non.cnv.events)
      } else {
        event.count <- nrow(events.explained)
      }
      
      total.frac <- event.count/total.events
      
      sample.cover[[i]] <- list(coverage.df = covered.sample.genomics, events.explained = events.explained,
                                event.count = event.count, total.frac = total.frac)
      
    }
    
    sample.cover
    
  })
  
  names(coverage) <- all.sample.names
  coverage
  
}

#' Merge Genomic Saturation 
#'  
#' Updated function to generate summary table of upstream genomics calculations
#' for plotting and determination of checkpoint
#' @param coverage.range list object generated by sampleEventOverlap function with per sample coverage information
#' @param topN number of regulators to consider
#' @return dataframe with summary information of cluster wide event coverage
#' @keywords internal
genomicSaturationSummary <- function(coverage.range, topN) {
  
  df <- tibble::tibble(.rows = 0)
  
  for (i in seq_len(topN)) {
    
    # get total event count for each sample at that k regulator
    count <- vapply(coverage.range, function(x) x[[i]]$event.count, numeric(1))
    count <- na.omit(count)
    
    # get all the fractions 
    fraction <- vapply(coverage.range, function(x) x[[i]]$total.frac, numeric(1))
    fraction <- na.omit(fraction)
    
    # to get all the unique events merge together the tables for each separate sample then get distinct
    unique.events <- purrr::map(coverage.range, list(i, "events.explained")) %>%
      dplyr::bind_rows() %>%
      dplyr::distinct() %>%
      nrow()
    
    data <- c(mean(count), i, mean(fraction), unique.events)
    
    df <- rbind(df, data)
  }
  
  colnames(df) <- c("mean", "k", "fraction", "unique.events")
  
  df
  
}

#' Get Saturation Inflection Point
#' 
#' Function to calculate checkpoint based on derivative of coverage fraction values
#' @param fractions vector of values representing fraction of events covered by that number of top regulators
#' @param clus.id cluster being evaluated
#' @return number representative of regulators in the checkpoint
#' @keywords internal
getInflection <- function(fractions, clus.id) {
  
  res <- pracma::gradient(fractions, h1 = 1)
  
  # get spot where there are at least 2 consecutive 0s in step gradient change
  zero.count <- 0
  index <- 0
  while(zero.count < 2 & index < length(res)) {
    index <- index + 1
    if(res[index] == 0){
      zero.count <- zero.count + 1
    } else {
      zero.count <- 0
    }
  } 
  
  if(index == length(res)) {
    index <- NA
  }

  if(is.na(index)) {
    message("Cluster ", clus.id, " did not saturate. Re-run analysis using a larger topN value")
    best <- 0
  } else {
    # mark inflection point as first in set of 1
    best <- index - 1
  }
  
  best
  
}


#' Function to generate null model for the saturation curves
#' @param momaObj momaObj 
#' @param viper.samples samples to consider
#' @param tissue.cluster which cluster to look at
#' @param topN number of MRs to take into account 
#' @param cytoband.collapse T/F regarding whether or not the genomic saturation analysis was done with collapsing the CNVs to cytobands
#' @param interaction.p pvalue threshold to use for MR-event interactions
#' @param aqtl.p pvalue threshold to use for aQTL scores
#' @param new T/F based on whether or not using results from version 1 or 2 of the pipeline
#' @return null means object
#' @keywords internal
saturationNull <- function(momaObj, viper.samples, tissue.cluster, topN, cytoband.collapse = T, interaction.p = 0.05, aqtl.p = 1, new = F) {
  message("Generating null model for saturation plots...")
  
  null.means <- c()
  
  # split analysis based on v1 vs v2
  if(isTRUE(new)){
    
    ranks <- momaObj$ranks.byCluster[[tissue.cluster]] %>% 
      dplyr::select(.data$regulator, .data$int.mr.p) %>%
      tibble::deframe()
    interactions.df <- dplyr::left_join(momaObj$interactions.byCluster[[tissue.cluster]], momaObj$gene.loc.mapping,
                                        by = c("event" = "Entrez.IDs")) %>%
      dplyr::filter(int.p <= interaction.p) %>%
      dplyr::filter(aQTL <= aqtl.p)
    
    for(i in seq_len(10)) {
      
      # look at the bottom half of TFs
      null.tfs <- ranks[which(ranks > quantile(ranks, 0.5))]
      set.seed(i)
      pvals <- sample(as.numeric(null.tfs), length(null.tfs))
      names(pvals) <- names(null.tfs)
      permuted.pvals <- sort(pvals)
      
      coverage.range <- sampleEventOverlap(momaObj, viper.samples, names(permuted.pvals),
                                           interactions.df, cytoband.collapse, topN, verbose = F)
      
      null.df <- genomicSaturationSummary(coverage.range, topN)
      null.means <- cbind(null.means, null.df$fraction)
    }
    
  } else {
    
    ranks <- momaObj$ranks[["integrated"]]
    
    for(i in seq_len(10)) {
      
      # look at the bottom half of TFs
      null.tfs <- ranks[which(ranks > quantile(ranks, 0.5))]
      set.seed(i)
      pvals <- sample(as.numeric(null.tfs), length(null.tfs))
      names(pvals) <- names(null.tfs)
      permuted.pvals <- sort(pvals)
      
      coverage.range <- getCoverage(momaObj, names(permuted.pvals), 
                                    viper.samples, topN = topN, verbose = F)
      
      null.df <- mergeGenomicSaturation(coverage.range, topN = topN)
      null.means <- cbind(null.means, null.df$fraction)
      
    }
  }
  
  null.means
}