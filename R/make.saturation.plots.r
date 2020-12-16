#' Main function to generate the summary plots of the analysis
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#' @param momaObj : momaObj that has already run the saturationCalculation function
#' @param clustering.solution : clustering vector with sample names and cluster designations
#' @param important.genes : vector of gene names to prioritize when plotting. 
#' Can be general genes of interest, oncogenes, tumor supressors etc
#' @param fCNV : vector of confirmed functional CNVs if calculated. Will filter for only those CNVs
#' @param max.events : maximum number of events to plot for the oncoplots
#' @return object with both types of summary plot for each subtype
#' @examples 
#' 
#' \dontrun{
#' makeSaturationPlots(momaObj, max.events = 20)
#' }
#' @export
makeSaturationPlots <- function(momaObj, clustering.solution = NULL, 
                                important.genes = NULL, fCNV = NULL, 
                                max.events = 30) {
  
  # check for required components of previous saturation analysis and set as local variables
  if (length(momaObj$checkpoints) == 0) {
    stop("Genomic Saturation analysis has not been done. Run momaObj$saturationCalculation() before continuing.
         Quitting...")
  }
  
  checkpoints <- momaObj$checkpoints
  genomic.saturation <- momaObj$genomic.saturation
  
  # get clustering solution to use for calculations. should be provided or previously saved in the momaObj
  if (is.null(clustering.solution)) {
    if (length(momaObj$sample.clustering) == 0) {
      stop("No clustering solution provided. Provide one as an argument or save one
                to the momaObj. Quitting...")
    } else {
      clustering.solution <- momaObj$sample.clustering
    }
  }
  
  # make gene to cytoband location mapping
  band2gene <- momaObj$gene.loc.mapping$Gene.Symbol
  names(band2gene) <- momaObj$gene.loc.mapping$Cytoband
  
  
  
  ###########################
  # Get raw input mutation/cnv/fusion matrices  
  # from the momaObj and prep for oncoprint plots   
  ##########################
  
  # Mutations
  snpmat <- momaObj$mut
  snpmat[snpmat == 1] <- "mut"
  snpmat[snpmat == 0] <- NA
  rownames(snpmat) <- mapEntrez(rownames(snpmat))
  
  # Fusions if they exist
  fusions.mat <- NULL
  if(!is.na(momaObj$fusions[1,1])) {
    fusions.mat <- momaObj$fusions
    fusions.mat[fusions.mat == 1] <- "fus"
    fusions.mat[fusions.mat == 0] <- NA
  }
  
  # CNVs. Input should be GISITC scores so separate into high/low events
  # Filter to only keep fCNVs if provided/in momaObj
  # TODO potentially change this to only use hypotheses results from momaobj?
  cnv <- momaObj$cnv
  if(!is.null(fCNV) | length(momaObj$fCNVs) > 0){
    cnv <- cnv[na.omit(match(fCNV, rownames(cnv))),]
    if(nrow(cnv) == 0) {
      stop("CNV matrix after functional CNV filtering is empty. 
           Check that names are in same format or if too many fCNVs are provided.
           Quitting...")
    }
  }
  
  amps <- dels <- cnv
  
  amps[amps < 0.5] <- NA
  amps[amps >= 0.5] <- "amp"
  
  dels[dels > -0.5] <- NA
  dels[dels <= -0.5] <- "del"
  
  
  rownames(amps) <- rownames(dels) <- mapEntrez(rownames(cnv))
  
  ##########################
  # Restructure saturation data 
  # to prep for plotting
  ##########################
  
  
  # get subtype event tables
  subtype.tables <- getSubtypeEventTables(genomic.saturation, clustering.solution, checkpoints)
  
  # get summary table of unique events added in for each regulator 
  ## could be clarified/improved
  ## also potential improvement: look for inflection points of huge jumps of new unique events and highlight those regulators in particular?
  tissue.coverage.df <- mergeDataBySubtype(genomic.saturation, clustering.solution, 100)
  
  
  # initalize object to save plots in
  tmp.oncoplots <- list()
  tmp.curveplots <- list()
  
  if(!is.null(important.genes)) {
    message("List of important genes has been provided. Will prioritize plotting these events first")
  }
  
  for(k in seq_along(subtype.tables)) {
    
    ###########
    # Oncoprint Events Plots
    ###########
    samples.thisCluster <- names(clustering.solution[clustering.solution == k])
    message("")
    message("Number of samples in cluster ", k, ": ", length(samples.thisCluster))
    
    # subset genomic event matrices for this cluster
    snpmat.thisClus <- snpmat[,colnames(snpmat) %in% samples.thisCluster]
    #cnv.thisClus <- cnv[,colnames(cnv) %in% samples.thisCluster]
    amps.thisClus <- amps[,colnames(amps) %in% samples.thisCluster]
    dels.thisClus <- dels[,colnames(dels) %in% samples.thisCluster]
    fusions.thisClus <- NULL
    if (!is.null(fusions.mat)) {	
      fusions.thisClus <- fusions.mat[,colnames(fusions.mat) %in% samples.thisCluster]
      if(ncol(fusions.thisClus) == 0) {
        fusions.thisClus <- NULL
      }
    }
    
    p.oncoprint <- oncoprintPlot(summary.vec = subtype.tables[[k]], snpmat.thisClus, amps.thisClus, dels.thisClus, fusions.thisClus, 
                                 important.genes, band2gene, max.events, k)
    tmp.oncoplots[[k]] <- p.oncoprint
    
    #########
    # Saturation Curve Plots
    #########
    p.coverage <- genomicPlotSmall(tissue.coverage.df, fraction=0.85, tissue.cluster=k)
    tmp.curveplots[[k]] <- p.coverage
    
  }
  
  list("oncoprint.plots" = tmp.oncoplots, "curve.plots" = tmp.curveplots)
  
}










#' Helper function to get subtype specific events
#' @param saturation.data : genomic saturation object from MOMA. List indexed by 
#' cluster then sample then regulator with the number of events associated with 
#' each additional regulator
#' @param sample.clustering : clustering vector with sample names and 
#' cluster designations
#' @param checkpoints : from momaObj
#' @return a table that has counts of how many times a particular event 
#' happens in a cluster
#' @keywords internal
getSubtypeEventTables <- function(saturation.data, sample.clustering, 
                                  checkpoints) {
  
  subtype.coverage <- list()
  coverage.allSubtypes <- saturation.data
  clustering <- sample.clustering
  
  # aggregate statistics across all samples, but use different subtype 
  # specific solutions
  for (cluster.id in unique(clustering)) {
    # dataframe with sample, type, gene names for each 
    # amp/del events can be locations. 
    coverage <- coverage.allSubtypes[[cluster.id]]
    sample.set <- names(coverage)
    
    # get number of regulators decided for this checkpoint depending on threshold 
    mr.cutoff <- length(checkpoints[[cluster.id]])
    
    #if (isTRUE(checkpoint.spec)) {
    #mr.cutoff <- length(momaObj$checkpoints.clusterSpecific[[cluster.id]])
    #} else {
    # mr.cutoff <- length(momaObj$checkpoints[[cluster.id]])
    #}
    
    subtype.coverage[[cluster.id]] <- makeCoverageDf(coverage, cutoff=mr.cutoff)
    
  }
  names(subtype.coverage) <- unique(clustering)
  
  # make summary table that counts the number of events in the group
  subtype.tables <- lapply(names(subtype.coverage), function(cluster.id) {
    table(apply(subtype.coverage[[cluster.id]], 1, function(x) paste0(x[1], ':', x[2])))
  })
  subtype.tables
}


#' Helper function for making the coverage dataframe 
#' @param coverage.list : List indexed by sample name, 
#' contains mut/fus/amp/del interactions
#' @param cutoff : number of regulators to include
#' @return dataframe with each sample and which events are 
#' captured by the checkpoint mrs
#' @keywords internal
makeCoverageDf <- function(coverage.list, cutoff) {
  
  df <- c()
  for (name in names(coverage.list)) {
    # skip empty/missing sample data
    if (is.null(coverage.list[[name]])) { next }
    
    for (type in c("mut", "amp", "del", "fus")) {
      
      events.thisSample <- coverage.list[[name]][[cutoff]][[type]]
      if (is.null(events.thisSample)) {
        non.null.idx <- which(vapply(coverage.list[[name]], 
                                     function(x) !is.null(x), 
                                     FUN.VALUE = logical(1))) ## ??
        real.cutoff <- max(non.null.idx[non.null.idx < cutoff])
        events.thisSample <- coverage.list[[name]][[real.cutoff]][[type]]
      }
      
      CT <- length(events.thisSample)
      if (CT==0) { next }
      submat <- cbind(rep(name, CT), rep(type, CT), events.thisSample)
      df <- rbind(df, submat)
    }
  }
  
  df <- data.frame(event=df[,3], type=df[,2], sample=df[,1])
  df	
}


#' Create data frame from coverage data, including number of total events 
#' 'covered' and unique events
#' @param genomic.saturation : data from genomic saturation function
#' @param sample.clustering : clustering vector with sample names and 
#' cluster designations
#' @param topN : number of regulators to look through. default is 100
#' @return dataframe with coverage data for genomic events
#' @keywords internal
mergeDataBySubtype <- function(genomic.saturation, sample.clustering, 
                               topN = 100)  {
  
  # generate summary stats for each subtype	
  full.df <- c()
  for (subtype in unique(sample.clustering)) {	
    
    # unnecessary 
    # subtype.samples <- names(sample.clustering[sample.clustering==subtype])
    coverage.subtype <- genomic.saturation[[subtype]]
    df <- mergeData(coverage.subtype, topN)
    
    # append a column specifying the subtype	
    df$subtype <- rep(subtype, nrow(df))
    full.df <- rbind(full.df, df)
  }
  
  full.df
}


#' Helper function for mergeDataBySubtype
#' @param coverage.range : genomic saturation for a particular subtype
#' @param topN : max number of top regulators to search through
#' @return dataframe with coverage data for genomic events
#' @keywords internal
mergeData <- function(coverage.range, topN)  {
  
  data <- c()
  for (i in seq_len(topN)) {
    # count for each sample
    # $mut/amp/del all point to either a NA or a vector of names of the event. If NA the length will be zero
    # so simply count the number of each type of event 
    count <- unlist(lapply(coverage.range, function(x) {
      num.events <- length(x[[i]]$mut)+length(x[[i]]$amp)+length(x[[i]]$del)
    }))
    count <- na.omit(count)
    
    # apply over each sample, get the coverage for each
    fraction <- unlist(lapply(coverage.range, function(x) {
      # critically: must omit the NAs so they don't interfere with count
      event.fractions <- x[[i]]$total.frac
      event.fractions
    }))
    fraction <- na.omit(fraction)
    
    all.events <- unlist(lapply(coverage.range, function(x) {
      c(x[[i]]$mut, x[[i]]$amp, x[[i]]$del)
    }))
    all.events <- na.omit(all.events)
    
    data <- rbind(data, c(i, mean(count), mean(fraction), length(unique(all.events))))
  }
  df <- data.frame(mean=data[,2], k=data[,1], fraction=data[,3], unique.events=data[,4]) 
  df	
}


#' Function to plot genomic events in the style of oncoPrint/cBioPortal
#' @importFrom tidyr separate
#' @importFrom dplyr arrange select mutate mutate_all everything bind_rows group_by n desc
#' @importFrom tibble deframe as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom stringr str_split_fixed
#' @importFrom rlang .data
#' @importFrom grid textGrob
#' @param summary.vec : named vector of the counts, named 'Event name':'Type'
#' where type is 'mut', 'amp', 'del', 'fus'. Mutations are in Entrez ID
#' Amp/Deletion CNV events are in genomic band location
#' @param snpmat.thisClus : SNP matrix subset to samples in current cluster
#' @param amps.thisClus : CNV matrix subset to samples in current cluster (just amplifications)
#' @param dels.thisClus : CNV matrix subset to samples in current cluster (just deletions)
#' @param fusions.thisClus : Fusion matrix subset to samples in current cluster
#' @param important.genes : well known genes to highlight in the analysis
#' @param band2gene : mapping of genomic location IDs to gene name: vector of 
#' HUGO gene ids, named by genomic location
#' @param max.events : maximum number of events to plot for the oncoplots
#' @param k : current cluster number
#' @return oncoprint event plot
#' @keywords internal
oncoprintPlot <- function(summary.vec, snpmat.thisClus, amps.thisClus, 
                          dels.thisClus, fusions.thisClus, important.genes, 
                          band2gene, max.events, k) {
  
  data <- data.frame(coverage=names(summary.vec), freq=as.numeric(summary.vec)) %>%
    tidyr::separate(.data$coverage, c("id", "type"), ":", remove = FALSE) %>%
    dplyr::arrange(dplyr::desc(.data$freq))
  
  ## Split into separate data types
  mut.data <- data[data$type == "mut",] %>% 
    dplyr::mutate(id = mapEntrez(.data$id)) %>% 
    dplyr::select(.data$id, .data$freq) %>% 
    tibble::deframe()
  
  amp.data <- data[data$type == "amp",] %>%
    dplyr::select(.data$id, .data$freq) %>% 
    dplyr::mutate(gene = NA)
  
  del.data <- data[data$type == "del",] %>%
    dplyr::select(.data$id, .data$freq) %>% 
    dplyr::mutate(gene = NA)
  
  fus.data <- data[data$type == "fus",] %>% 
    dplyr::select(.data$id, .data$freq) %>% 
    tibble::deframe()
  
  
  if(is.null(important.genes)) {
    # if no list of important genes has been provided just plot top events
    
    # subset mutations by name
    mut.mat <- tibble::as_tibble(snpmat.thisClus[names(mut.data),], rownames = NA) %>% 
      tibble::rownames_to_column(var = "genomic.event") %>% 
      dplyr::mutate(type = "mut") %>% 
      dplyr::mutate_all(as.character) %>% 
      dplyr::select(.data$genomic.event, .data$type, dplyr::everything())
    
    # for each cytoband region subset the genes in that region and get the gene 
    # with the highest number of events to represent it
    
    # amplifications 
    for(row in seq_len(nrow(amp.data))) {
      cnv.loc <- amp.data$id[row]
      genes.inBand <- as.character(band2gene[which(names(band2gene)==cnv.loc)])
      band.mat <- amps.thisClus[intersect(genes.inBand, rownames(amps.thisClus)),]
      if(is.null(nrow(band.mat))) { 
        # only one gene case, just use that one
        amp.data$gene[row] <- intersect(genes.inBand, rownames(amps.thisClus))
      } else {
        rank.events <- apply(band.mat, 1, function(x) sum(!is.na(x))) %>% sort(decreasing = TRUE)  
        # take the one with most events cause sometimes genes in the band won't all have same score
        amp.data$gene[row] <- names(rank.events)[1]
      }
      
    }
    
    amps.mat <- amps.thisClus[amp.data$gene, ]
    amps.mat <- cbind(genomic.event = amp.data$id, type = rep("amp", nrow(amps.mat)), amps.mat) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    # deletions
    for(row in seq_len(nrow(del.data))) {
      cnv.loc <- del.data$id[row]
      genes.inBand <- as.character(band2gene[which(names(band2gene)==cnv.loc)])
      band.mat <- dels.thisClus[intersect(genes.inBand, rownames(dels.thisClus)),]
      if(is.null(nrow(band.mat))) { 
        # only one gene case, just use that one
        del.data$gene[row] <- intersect(genes.inBand, rownames(dels.thisClus))
      } else {
        rank.events <- apply(band.mat, 1, function(x) sum(!is.na(x))) %>% sort(decreasing = TRUE)  
        # take the one with most events cause sometimes genes in the band won't all have same score
        del.data$gene[row] <- names(rank.events)[1]
      }
      
    }
    
    dels.mat <- dels.thisClus[del.data$gene, ]
    dels.mat <- cbind(genomic.event = del.data$id, type = rep("del", nrow(dels.mat)), dels.mat) %>%
      as.data.frame(stringsAsFactors = FALSE)
    
    # fusions
    
    if(!is.null(fusions.thisClus) & length(fus.data) > 0) {
      fus.mat <- tibble::as_tibble(fusions.thisClus[names(fus.data),], rownames = NA) %>% 
        tibble::rownames_to_column(var = "genomic.event") %>% 
        dplyr::mutate(type = "fus") %>% 
        dplyr::mutate_all(as.character) %>% 
        dplyr::select(.data$genomic.event, .data$type, dplyr::everything())
    } else {
      fus.mat <- NULL
    }
    
    
    mat.to.plot <- dplyr::bind_rows(mut.mat, amps.mat, dels.mat, fus.mat)
    
  } 
  
  
  # merge rows that are duplicated (have muts and amp/dels)
  dups.count <- mat.to.plot %>% dplyr::group_by(.data$genomic.event) %>% 
    dplyr::summarise(n = dplyr::n())
  dups.count <- dups.count[dups.count$n > 1, 1]
  dups <- mat.to.plot[mat.to.plot$genomic.event %in% dups.count$genomic.event,]
  
  
  if(nrow(dups) > 0 ) {
    dups.names <- unique(dups$genomic.event)
    # remove duplicates from the current matrix
    mat.to.plot <- mat.to.plot[!mat.to.plot$genomic.event %in% dups.count$genomic.event,]
    for (name in dups.names) {
      # collapse all types of events together into one row and then rebind to main matrix
      # having NA;event is fine syntax for use in oncoprint function
      #### added in function to keep amps and dels separate, only merge mut with amps/dels
      
      merged <- dups[dups$genomic.event == name,]
      if(!"mut" %in% merged$type) {
        # means that only events are amps/dels, don't merge. 
        # add a/d to the end do be able to differentiate later 
        # (and to not run into issue of duplicate rownames in the matrix)
        merged$genomic.event <- paste(merged$genomic.event, merged$type, sep = "_")
      } else if ("mut" %in% merged$type & nrow(merged) == 2) {
        # means that one event is a mut and other is amp/del
        # merge together like before by collapsing the rows together
        merged <- apply(merged, 2, paste0, collapse = ";")
        merged[1] <- name
      } else if (nrow(merged) == 3) {
        # means that there is one of each type of event mut, amp and del
        # merge mut with amp and del separately
        
        amp.merged <- merged[merged$type %in% c("amp", "mut"),] %>% apply(2, paste0, collapse = ";")
        amp.merged[1] <- paste(name, "amp", sep = "_")
        
        del.merged <- merged[merged$type %in% c("del", "mut"),] %>% apply(2, paste0, collapse = ";")
        del.merged[1] <- paste(name, "del", sep = "_")
        
        merged <- dplyr::bind_rows(amp.merged, del.merged)
      } else {
        stop("Problem merging duplicates, double check gene names.")
      }
      
      mat.to.plot <- dplyr::bind_rows(mat.to.plot, merged)
    }
  }
  
  # replace new NA's with regular ones so that plots can be counted for events again
  mat.to.plot[mat.to.plot == "NA;NA"] <- NA
  mat.to.plot[mat.to.plot == "NA;NA;NA"] <- NA
  
  # get row sums, in terms of events and take top ones. filter for events that happen at least twice here
  mat.to.plot <- mat.to.plot %>% dplyr::select(-.data$type) %>%
    tibble::column_to_rownames(var = "genomic.event") %>% 
    as.matrix()
  num.events <- apply(mat.to.plot, 1, function(x) sum(!is.na(x))) %>% sort(decreasing = TRUE) 
  num.events <- num.events[num.events > 1] 
  
  # add events in chunks of frequency until max is reached
  final.mat <- matrix(nrow = 0, ncol = ncol(mat.to.plot))
  freq.range <- unique(num.events)
  freq.idx <- 1
  while(nrow(final.mat) < max.events & freq.idx <= length(unique(num.events))) {
    to.add <- names(num.events[num.events == freq.range[freq.idx]])
    to.add.mat <- matrix(mat.to.plot[to.add,], nrow = length(to.add), dimnames = list(rownames = to.add))
    # final.mat <- rbind(final.mat, to.add.mat)
    new.final.mat <- rbind(final.mat, to.add.mat)
    if(nrow(new.final.mat) > max.events) {
      break
    } else {
      final.mat <- new.final.mat
    }
    freq.idx <- freq.idx + 1
  }
  
  # check if all genes are unique (as per merging of muts with amps/dels)
  # if so remove _cnv ending, else leave it on
  
  final.plot.names <- stringr::str_split_fixed(rownames(final.mat), pattern = "_", n = 2)
  if(length(unique(final.plot.names[,1])) == nrow(final.mat)) {
    rownames(final.mat) <- final.plot.names[,1]
  }
  
  message("Plotting ", nrow(final.mat), " total events for subtype ")
  
  if(nrow(final.mat) == 0) {
    return(grid::textGrob("No significant events to plot for this subtype"))
  }
  
  
  ##########
  # Heatmap plot
  #########
  
  
  heatmap_legend_param = list(title = "Alterations", 
                              at = c("amp", "del", "mut", "fus"), 
                              labels = c("Amplification", "Deletion", "Mutation", "Fusion"))
  
  
  # color parameters for oncoprint
  col = c("del" = "blue", "amp" = "red", "mut" = "#008000", "fus" = 'gold2', "NA" = "grey22")
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w*0.6, h*0.75, 
                gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # big blue - dels
    del = function(x, y, w, h) {
      grid.rect(x, y, w*0.6, h*0.75, 
                gp = gpar(fill = col["del"], col = NA))
    },
    # big red - amps
    amp = function(x, y, w, h) {
      grid.rect(x, y, w*0.6, h*0.75, 
                gp = gpar(fill = col["amp"], col = NA))
    },
    
    # small green
    mut = function(x, y, w, h) {
      grid.rect(x, y, w*0.6, h*0.33, 
                gp = gpar(fill = col["mut"], col = NA))
    },
    
    # big yellow
    fus = function(x, y, w, h) {
      grid.rect(x, y, w*0.6, h*0.75, 
                gp = gpar(fill = col["fus"], col = NA))
    }
  )
  
  # optimize size of labels
  if(nrow(final.mat) >= 55) {
    label.size <- 3
  } else if (nrow(final.mat) >= 40 ) {
    label.size <- 4
  } else if (nrow(final.mat) >= 30) {
    label.size <- 5
  } else if (nrow(final.mat) >= 20) {
    label.size <- 6
  } else {
    label.size <- 7
  }
  
  
  title <- paste("Genomic Events for Subtype", k)
  
  # set NAs to empty string in final mat to avoid Complex Heatmap warning
  final.mat[is.na(final.mat)] <- ""
  
  ht <- ComplexHeatmap::oncoPrint(final.mat,
                                  alter_fun = alter_fun, col = col,
                                  heatmap_legend_param = heatmap_legend_param,
                                  pct_side = "right", pct_gp = gpar(fontsize = label.size),
                                  column_title = title, column_title_gp = gpar(fontsize = 8),
                                  row_names_side = "left", row_names_gp = gpar(fontsize = label.size),
                                  show_heatmap_legend = TRUE,
                                  top_annotation = HeatmapAnnotation(
                                    column_barplot = anno_oncoprint_barplot(height = unit(.4, "cm"), axis_param = list(gp = gpar(fontsize = 4)))),
                                  right_annotation = rowAnnotation(
                                    row_barplot = anno_oncoprint_barplot(width = unit(1, "cm"), axis_param = list(gp = gpar(fontsize = 5)))))
  
  
  
  p <- grid::grid.grabExpr(draw(ht, padding = unit(c(0, 0, 0, 0), "mm")))
  
  p
  
}







##### No longer the primary plot style but keeping as an option plot style
##### TODO: Need to make it a freestanding function
#' Plot barchart of genomic events 
#' 
#' @importFrom rlang .data
#' @param summary.vec : named vector of the counts, named 'Event name':'Type'
#' where type is 'mut', 'amp', 'del', 'fus'. Mutations are in Entrez ID
#' Amp/Deletion CNV events are in genomic band location
#' @param highlight.genes : well known genes to highlight in the analysis in 
#' @param genomeBand_2_gene : mapping of genomic location IDs to gene name: 
#' vector of HUGO gene ids, named by genomic loci
#' @param samples.total : number of samples in the subtype, used to calculate percentages
#' @param max.muts : maximum number of mutations to get per sample, default is 10
#' @param max.cnv : maximum number of cnvs to per sample, default is 5
#' @return plot object
#' @keywords internal 
plotEvents <- function(summary.vec, highlight.genes=NULL, genomeBand_2_gene=NULL,
                       samples.total, max.muts = 10, max.cnv = 5) {
  
  data <- data.frame(coverage=names(summary.vec), Freq=as.numeric(summary.vec))
  # order by individual frequency
  data <- data[order(-data$Freq),]
  
  data$id <- unlist(lapply(data$coverage, function(label) {
    label <- as.character(label)
    name <- strsplit(label, ':')[[1]][1]
    hugo <- mapEntrez(as.character(name))
    if (is.na(hugo)) {
      return (name)
    } else {
      return (hugo)
    }
  }))
  
  # make a new column: for CNV band locations
  data$type <- unlist(lapply(data$coverage, function(label) {
    type <- strsplit(as.character(label), ':')[[1]][2]
    type
  }))
  
  # get order dataframe of events 
  mapped <- getDataFrame(data, highlight.genes, genomeBand_2_gene, max.muts = 10, max.cnv = 5)
  # check the size: if we have too many events to display, then select the top N unique IDs
  # already sorted by frequency of occurence. 
  if (nrow(mapped) > 50) {
    # add at least half simply with mutated gene labels
    mut.data <- mapped[mapped$type=='mut',]
    add.labels <- na.omit(unique(mut.data[order(-mut.data$Freq),][seq_len(25),]$id))
    # and add cnv
    additional <- setdiff(unique(mapped$id), add.labels)[seq_len(50-length(add.labels))]
    mapped <- mapped[mapped$id %in% union(additional, add.labels),]
  }
  
  # The palette with black: (unnecessary?)
  # cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # if exporting as object then just use default of 10 and then change by adding geom/layer when plotting
  # else if plotting directly then try to adjust text size to fit based on number of events
  
  # data.size <- nrow(mapped)
  # if(is.null(output)) {
  #   print("No file output name found, using pre-selected text size for plotting")
  #   y.textSize <- 10
  # } else {
  #   print("Output plot name found, adjusting text size based on number of events")
  #   if (num.subtypes > 5 && data.size > 50) {
  #     y.textSize <- 1
  #   } else if (num.subtypes > 4 && data.size > 25) {
  #     y.textSize <- 2
  #   } else if (num.subtypes > 4 && data.size <= 25) {
  #     y.textSize <- 3
  #   } else if (data.size < 20) {
  #     y.textSize <- 5
  #   } else if (data.size < 40) {
  #     y.textSize <- 6
  #   }
  # }
  
  message('Number of entries in events matrix: ', nrow(mapped))
  #print (paste('Using font size ', y.textSize))
  
  
  # Scale frequency to a percentage of samples in that cluster not number of occurences
  mapped$freq.percentage <- mapped$Freq/samples.total
  
  p <- ggplot2::qplot(x=mapped$id, y=mapped$freq.percentage, fill=mapped$type, data=mapped, geom = "col") + coord_flip() +
    scale_fill_manual(values = c("mut"='#00BA38', "amp"= '#F8766D', "del" = '#619CFF', "fus" = '#FF8C00' )) +
    ylab("Frequency") + xlab("Event") 
  #theme(axis.text.y = element_text(size=y.textSize))
  #labs(title=label) 
  
  # if (!is.null(output)) {
  #   p <- p + labs(title = paste(plot.label, "subtype", k))
  #   ggsave(output)
  # }
  
  p
}



#' Helper function to get data frame for bar plot plot.events function

#' @param data data.frame with $type, $id, $Freq per event
#' @param highlight.genes genes to look for in mutations/cnv lists (if looking
#'  for specific genes because of prior knowledge)
#' @param genomeBand_2_gene mapping of genomic location IDs to gene name: 
#' vector of HUGO gene ids, named by genomic loci
#' @param max.muts maximum number of mutations to get per sample, default is 10
#' @param max.cnv maximum number of cnvs to per sample, default is 5
#' @return ordered data frame with each genomic event and it's frequency
#' @keywords internal
getDataFrame <- function(data, highlight.genes, genomeBand_2_gene, 
                         max.muts = 10, max.cnv = 5) {
  
  #### edit logic of this section ?
  
  
  MAX.muts <- max.muts
  MAX.CNV.loc <- max.cnv
  ## all all highlight genes and events to the 
  
  ## for each CNV location event, find genes that overlap with the highlight list. 
  ## add duplicate entries for those genes
  mapped <- c()
  loc.data <- data[data$type=='del',]
  for (row in seq_len(nrow(loc.data))) {
    loc <- loc.data[row, 3]
    genes.inBand <- as.character(genomeBand_2_gene[which(names(genomeBand_2_gene)==loc)])
    hgIB <- intersect(genes.inBand, highlight.genes)
    
    if (length(hgIB)==0) { next }
    
    mapped <- rbind(mapped, data.frame(coverage=loc.data[row, 1], 
                                       Freq=loc.data[row, 2], 
                                       id=hgIB, 
                                       type=loc.data[row, 4]))
  }
  loc.data <- data[data$type=='amp',]
  for (row in seq_len(nrow(loc.data))) {
    loc <- loc.data[row, 3]
    genes.inBand <- as.character(genomeBand_2_gene[which(names(genomeBand_2_gene)==loc)])
    hgIB <- intersect(genes.inBand, highlight.genes)
    
    if (length(hgIB)==0) { next }
    
    mapped <- rbind(mapped, data.frame(coverage=loc.data[row, 1], 
                                       Freq=loc.data[row, 2], 
                                       id=hgIB, 
                                       type=loc.data[row, 4]))
  }
  
  # find mutations in key regions
  loc.data <- data[data$type=='mut',]
  for (row in seq_len(nrow(loc.data))) {
    gene <- loc.data[row, 3]
    hgIB <- intersect(gene, highlight.genes)
    if (length(hgIB)==0) { next }
    mapped <- rbind(mapped, data.frame(coverage=loc.data[row, 1], 
                                       Freq=loc.data[row, 2], 
                                       id=hgIB, 
                                       ype=loc.data[row, 4]))
  }
  
  # HUGO ids: add additional mutation events outside of the highlight genes
  all.gene.ids <- as.character(mapped$id)
  # add in top X mutations not in the driver list
  mut.data <- data[data$type=='mut',]
  i = 1	
  for (row in seq_len(nrow(mut.data))) {
    if (mut.data[row,]$id %in% all.gene.ids) { next }
    if (i > MAX.muts) { break }
    mapped <- rbind(mapped, mut.data[row,])
    i = 1+i
  }
  
  # add all fusions
  mapped <- rbind(mapped, data[data$type=='fus',])
  
  # add CNV band locations (a few)
  subset <- data[apply(cbind(data$type=='amp', data$type=='del'), 1, any),]
  if (nrow(subset) > MAX.CNV.loc) { subset <- subset[seq_len(MAX.CNV.loc),] }
  mapped <- rbind(mapped, subset)
  
  # order by total frequency of events by gene/id
  mapped$event_sums <- vapply(mapped$id, function(id) { 
    sum(as.numeric(mapped[which(mapped$id==as.character(id)),]$Freq)) }, 
    FUN.VALUE = numeric(1))
  mapped <- mapped[order(-mapped$event_sums),]
  mapped$id <- factor(mapped$id, levels=unique(rev(mapped$id)))
  
  mapped
}


#' Make small genomic plot
#' @importFrom tidyr drop_na
#' @importFrom rlang .data
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#' @import magrittr
#' @param input.df : tissue.coverage.df with mean, k, fraction and unique events.
#' @param fraction : what fraction coverage to use for genomic curve threshold
#' @param tissue.cluster : which cluster subsample to look at
#' @return output .png
#' @keywords internal
genomicPlotSmall <- function(input.df, fraction=0.85, tissue.cluster=NULL)  { 
  
  #### need to add in null matrix function ###
  
  
  # get number of subtypes and colors for plotting
  num.subtypes <- length(unique(input.df$subtype))
  getPalette <- colorRampPalette(brewer.pal(n = 8, name = "Dark2"))
  subtype.colors <- getPalette(num.subtypes)
  color.this.sub <- subtype.colors[tissue.cluster]
  
  
  subtype.df <- input.df[input.df$subtype == tissue.cluster,] 
  subtype.df <- subtype.df %>% tidyr::drop_na(.data$fraction)
  
  
  sweep <- subtype.df$fraction
  names(sweep) <- sort(unique(subtype.df$k))
  best.k <- fitCurvePercent(sweep, frac=fraction)
  message("Number of MRs in checkpoint: ", best.k)
  
  
  # mean statistic for these samples
  mean.stat <- unlist(lapply(sort(unique(subtype.df$k)), function(k) {
    mean(subtype.df[which(subtype.df$k==k),]$fraction)
  }))
  sweep <- mean.stat
  # Get mean # of samples
  mean.events.stat <- unlist(lapply(sort(unique(subtype.df$k)), function(k) {
    mean(subtype.df[which(subtype.df$k==k),]$mean)
  }))
  names(mean.events.stat) <- sort(unique(input.df$k))
  # the average multiplier across all the bands:
  # use to count the mean number of events on the right side
  y.multiplier <- mean(na.omit(mean.events.stat/mean.stat))
  
  # first make a condensed plot just the first ~100 events
  ymax <- subtype.df[nrow(subtype.df),]$fraction
  
  # make plot
  p.100 <- ggplot(subtype.df, aes(.data$k, .data$fraction)) + geom_line(color = color.this.sub, size=1.5, alpha=0.75) +
    xlab("Number of MRs") +
    scale_y_continuous(
      "Mean Fraction",
      sec.axis = sec_axis(~ . * y.multiplier, name = "Count"),
      limits=c(0,ymax+0.01)
    ) +
    #geom_ribbon(aes(ymin=null.bq, ymax=null.uq), color='#808080', alpha=0.2, linetype=2) +
    geom_vline(xintercept=as.numeric(best.k), linetype=3) +
    xlim(0,100) +
    #labs(title=tissue) +
    theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), legend.position="none")
  
  # return the plot itself
  p.100 
}


