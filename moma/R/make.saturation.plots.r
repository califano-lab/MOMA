#' Helper function to get subtype specific events
#' @param saturation.data : genomic saturation object from moma. List indexed by cluster then sample then regulator with the number of events associated with each additional regulator
#' @param sample.clustering : clustering vector with sample names and cluster designations
#' @param checkpoints : from momaObj
#' @return a table that has counts of how many times a particular event happens in a cluster
get.subtype.event.tables <- function(saturation.data, sample.clustering, checkpoints) {
  
  subtype.coverage <- list()
  coverage.allSubtypes <- saturation.data
  clustering <- sample.clustering
  
  # aggregate statistics across all samples, but use different subtype specific solutions
  for (cluster.id in unique(clustering)) {
    # dataframe with sample, type, gene names for each 
    # amp/del events can be locations. 
    coverage <- coverage.allSubtypes[[cluster.id]]
    sample.set <- names(coverage)
    
    # get number of regulators decided for this checkpoint depending on coverage threshold 
    mr.cutoff <- length(checkpoints[[cluster.id]])
    
    #if (isTRUE(checkpoint.spec)) {
      #mr.cutoff <- length(momaObj$checkpoints.clusterSpecific[[cluster.id]])
    #} else {
     # mr.cutoff <- length(momaObj$checkpoints[[cluster.id]])
    #}
    
    subtype.coverage[[cluster.id]] <- make.coverage.df(coverage, cutoff=mr.cutoff)
    
  }
  names(subtype.coverage) <- unique(clustering)
  
  # make summary table that counts the number of events in the group
  subtype.tables <- lapply(names(subtype.coverage), function(cluster.id) {
    table(apply(subtype.coverage[[cluster.id]], 1, function(x) paste0(x[1], ':', x[2])))
  })
  subtype.tables
}


#' Helper function for making the coverage dataframe (might)
#' @param coverage.list : List indexed by sample name, contains $mut, $fus, $amp, $del interactions
#' @param cutoff : number of regulators to include
#' @return dataframe with each sample and which events are captured by the checkpoint mrs
make.coverage.df <- function(coverage.list, cutoff) {

df <- c()
for (name in names(coverage.list)) {
  # skip empty/missing sample data
  if (is.null(coverage.list[[name]])) { next }
  
  for (type in c("mut", "amp", "del", "fus")) {
    
    events.thisSample <- coverage.list[[name]][[cutoff]][[type]]
    if (is.null(events.thisSample)) {
      non.null.idx <- which(sapply(coverage.list[[name]], function(x) !is.null(x))) ### not sure that this is for
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


#' Create data frame from coverage data, including number of total events 'covered' and unique events
#' @param genomic.saturation : data from genomic saturation function
#' @param sample.clustering : clustering vector with sample names and cluster designations
#' @param topN : number of regulators to look through. default is 100
#' @return
merge.data.by.subtype <- function(genomic.saturation, sample.clustering, topN = 100)  {
  
  ### unnecessary ### remove after testing
  # flatten the coverage into a single list for all samples
  # coverage.allSamples <- list()
  # for (clus in 1:length(genomic.saturation)) {
  #   ll <- genomic.saturation[[clus]]
  #   for (sample in names(ll)) {
  #     coverage.allSamples[[sample]] <- ll[[sample]]
  #   }
  # }
  
  # generate summary stats for each subtype	
  full.df <- c()
  for (subtype in unique(sample.clustering)) {	
    
    # unnecessary 
    # subtype.samples <- names(sample.clustering[sample.clustering==subtype])
    coverage.subtype <- genomic.saturation[[subtype]]
    df <- merge.data(coverage.subtype, topN)
    
    # append a column specifying the subtype	
    df$subtype <- rep(subtype, nrow(df))
    full.df <- rbind(full.df, df)
  }
  
  full.df
}


#' Helper function for merge.data.by.subtype
#' @param coverage.range : genomic saturation for a particular subtype
#' @param topN : max number of top regulators to search through
#' @return
merge.data <- function(coverage.range, topN)  {
  
  data <- c()
  for (i in 1:topN) {
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


#' Plot barchart of genomic events
#'
#' @param summary.vec : named vector of the counts, named 'Event name':'Type'
#' where type is 'mut', 'amp', 'del', 'fus'. Mutations are in Entrez ID
#' Amp/Deletion CNV events are in genomic band location
#' @param highlight.genes : well known genes to highlight in the analysis in 
#' @param genomeBand_2_gene : mapping of genomic location IDs to gene name: vector of HUGO gene ids, named by genomic loci
#' @param samples.total : number of samples in the subtype, used to calculate percentages
#' @param max.muts : maximum number of mutations to get per sample, default is 10
#' @param max.cnv : maximum number of cnvs to per sample, default is 5
#' @return plot object 
plot.events <- function(summary.vec, highlight.genes=NULL, genomeBand_2_gene=NULL, samples.total, max.muts = 10, max.cnv = 5) {
  
  data <- data.frame(coverage=names(summary.vec), Freq=as.numeric(summary.vec))
  # order by individual frequency
  data <- data[order(-data$Freq),]
  
  data$id <- unlist(lapply(data$coverage, function(label) {
    label <- as.character(label)
    name <- strsplit(label, ':')[[1]][1]
    hugo <- map.entrez(as.character(name))
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
  mapped <- get.data.frame(data, highlight.genes, genomeBand_2_gene, max.muts = 10, max.cnv = 5)
  # check the size: if we have too many events to display, then select the top N unique IDs
  # already sorted by frequency of occurence. 
  if (nrow(mapped) > 50) {
    # add at least half simply with mutated gene labels
    mut.data <- mapped[mapped$type=='mut',]
    add.labels <- na.omit(unique(mut.data[order(-mut.data$Freq),][1:25,]$id))
    # and add cnv
    additional <- setdiff(unique(mapped$id), add.labels)[1:(50-length(add.labels))]
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
  
  print (paste('Number of entries in events matrix: ', nrow(mapped)))
  #print (paste('Using font size ', y.textSize))
  
  
  # Scale frequency to a percentage of samples in that cluster not number of occurences
  mapped$freq.percentage <- mapped$Freq/samples.total
  
  p <- ggplot2::qplot(x=id, y=freq.percentage, fill=type, data=mapped, geom = "col") + coord_flip() +
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

#' @param data : data.frame with $type, $id, $Freq per event
#' @param highlight.genes : genes to look for in mutations/cnv lists (if looking for specific genes because of prior knowledge)
#' @param genomeBand_2_gene : mapping of genomic location IDs to gene name: vector of HUGO gene ids, named by genomic loci
#' @param max.muts : maximum number of mutations to get per sample, default is 10
#' @param max.cnv : maximum number of cnvs to per sample, default is 5
#' @return ordered data frame with each genomic event and it's frequency
get.data.frame <- function(data, highlight.genes, genomeBand_2_gene, max.muts = 10, max.cnv = 5) {
  
  #### edit logic of this section ?
  
  
  MAX.muts <- max.muts
  MAX.CNV.loc <- max.cnv
  ## all all highlight genes and events to the 
  
  ## for each CNV location event, find genes in that overlap with the highlight list. 
  ## add duplicate entries for those genes
  mapped <- c()
  loc.data <- data[data$type=='del',]
  for (row in 1:nrow(loc.data)) {
    loc <- loc.data[row, 3]
    genes.inBand <- as.character(genomeBand_2_gene[which(names(genomeBand_2_gene)==loc)])
    hgIB <- intersect(genes.inBand, highlight.genes)
    
    if (length(hgIB)==0) { next }
    
    mapped <- rbind(mapped, data.frame(coverage=loc.data[row, 1], Freq=loc.data[row, 2], id=hgIB, type=loc.data[row, 4]))
  }
  loc.data <- data[data$type=='amp',]
  for (row in 1:nrow(loc.data)) {
    loc <- loc.data[row, 3]
    genes.inBand <- as.character(genomeBand_2_gene[which(names(genomeBand_2_gene)==loc)])
    hgIB <- intersect(genes.inBand, highlight.genes)
    
    if (length(hgIB)==0) { next }
    
    mapped <- rbind(mapped, data.frame(coverage=loc.data[row, 1], Freq=loc.data[row, 2], id=hgIB, type=loc.data[row, 4]))
  }
  
  # find mutations in key regions
  loc.data <- data[data$type=='mut',]
  for (row in 1:nrow(loc.data)) {
    gene <- loc.data[row, 3]
    hgIB <- intersect(gene, highlight.genes)
    if (length(hgIB)==0) { next }
    mapped <- rbind(mapped, data.frame(coverage=loc.data[row, 1], Freq=loc.data[row, 2], id=hgIB, type=loc.data[row, 4]))
  }
  
  # HUGO ids: add additional mutation events outside of the highlight genes
  all.gene.ids <- as.character(mapped$id)
  # add in top X mutations not in the driver list
  mut.data <- data[data$type=='mut',]
  i = 1	
  for (row in 1:nrow(mut.data)) {
    if (mut.data[row,]$id %in% all.gene.ids) { next }
    if (i > MAX.muts) { break }
    mapped <- rbind(mapped, mut.data[row,])
    i = 1+i
  }
  
  # add all fusions
  mapped <- rbind(mapped, data[data$type=='fus',])
  
  # add CNV band locations (a few)
  subset <- data[apply(cbind(data$type=='amp', data$type=='del'), 1, any),]
  if (nrow(subset) > MAX.CNV.loc) { subset <- subset[1:MAX.CNV.loc,] }
  mapped <- rbind(mapped, subset)
  
  # order by total frequency of events by gene/id
  mapped$event_sums <- sapply(mapped$id, function(id) { sum(as.numeric(mapped[which(mapped$id==as.character(id)),]$Freq)) })
  mapped <- mapped[order(-mapped$event_sums),]
  mapped$id <- factor(mapped$id, levels=unique(rev(mapped$id)))
  
  mapped
}


#' Make small genomic plot
#' @importFrom dplyr filter
#' @importFrom tidyr drop_na
#' @import ggplot2
#' @param input.df : tissue.coverage.df with mean, k, fraction and unique events. has all samples
#' @param fraction : what fraction coverage to use for genomic curve threshold
#' @param tissue.cluster : which cluster subsample to look at
#' @return output .png
genomic.plot.small <- function(input.df, fraction=0.85, tissue.cluster=NULL)  { 
  
  #### need to add in null matrix function? ###
  
  # old color code
  # subtype.color.map=brewer.pal(n = 8, name = "Dark2")
  # names(subtype.color.map) <- 1:8
  
  
  # get number of subtypes and colors for plotting
  num.subtypes <- length(unique(input.df$subtype))
  getPalette <- colorRampPalette(brewer.pal(n = 8, name = "Dark2"))
  subtype.colors <- getPalette(num.subtypes)
  color.this.sub <- subtype.colors[tissue.cluster]
  
  # subset to only this subtype and remove NAs
  subtype.df <- input.df %>% 
    dplyr::filter(subtype == tissue.cluster) %>%
    tidyr::drop_na(fraction)
  
  # df <- data.frame(k=input.df$k, mean=input.df$fraction, subtype=as.factor(input.df$subtype))
  
  sweep <- subtype.df$fraction
  names(sweep) <- sort(unique(subtype.df$k))
  best.k <- fit.threshold(sweep, frac=fraction)
  print (paste("Threshold for MR cutoff: ", best.k))
  
  
  ### this needs to be cleaned up too.... ####
  
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
  
  # need to fix so that colors are different for each plot
  # for now it's just black
  p.100 <- ggplot(subtype.df, aes(k, fraction)) + geom_line(color = color.this.sub, size=1.5, alpha=0.75) +
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


#'  Helper function to do inflection point fitting. 
#'  Use the heuristic, or if supplied, the definitions file to set the 
#'  MR threshold
#'  @param sweep : numeric vector of genomic coverage values, named by -k- threshold
#'  @param frac : Fraction of coverage to use as a threshold (default .85 = 85 percent)
#'  @return The -k- integer where coverage is acheived
fit.threshold <- function(sweep, frac=NULL) {
  k <- fit.curve.percent(sweep, frac)
  return (k)
}


# # Make multiplot layout
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   library(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }

