
#' @title MOMA Runner
#' @description Main class encapsulating the input data and logic of the MOMA algorithm
#' @import stats
#' @import qvalue
#' @import parallel
#' @import reshape2
#' @import MKmisc
#' @import survival
#' @import RColorBrewer
#' @import methods
#' @import ggplot2
#' @import magrittr
#' @import graphics
#' @import grDevices
#' @import readr
#' @import tibble
#' @export
momaRunner <- setRefClass("momaRunner", fields = 
                            list(viper = "matrix", 
                                 mut = "matrix", 
                                 cnv = "matrix", 
                                 fusions = "matrix", 
                                 pathways = "list", 
                                 gene.blacklist = "character",
                                 output.folder = "character", 
                                 gene.loc.mapping = "data.frame", 
                                 nes = "list", # result field
                                 interactions = "list", # result field
                                 clustering.results = "list", 
                                 ranks = "list", 
                                 hypotheses = "list",
                                 genomic.saturation = "list", 
                                 coverage.summaryStats = "list", 
                                 checkpoints = "list", 
                                 sample.clustering = "numeric", # numbers are cluster assignments, names are sample ids matching other data
                                 identity.plots = "list"), # result field
                          methods = list(
  runDIGGIT = function(fCNV = NULL, cnvthr = 0.5, min.events = 4, verbose = FALSE) {
      
      cnv.local <- NULL
      if (is.null(fCNV)) {
          print("Warning: no fCNV supplied, using no CNV filter!")
          cnv.local <- cnv
      } else {
          cnv.local <- cnv[intersect(fCNV, rownames(cnv)), ]
      }
      
      somut <- mut
      
      # Define amplifications and deletions
      amps <- dels <- cnv.local
      amps[which(amps < cnvthr)] <- 0
      amps[which(amps >= cnvthr)] <- 1
      dels[which(dels > -cnvthr)] <- 0
      dels[which(dels <= -cnvthr)] <- 1
      
      # save the exact hypotheses (genes) we're testing, based on MIN.EVENTS
      
      temp.amps <- rownames(amps[apply(amps, 1, sum, na.rm = TRUE) >= min.events, ])
      amps.hypotheses <- temp.amps[which(!(temp.amps %in% gene.blacklist))]
      amps.mat <- amps[amps.hypotheses,]
      
      temp.dels <- rownames(dels[apply(dels, 1, sum, na.rm = TRUE) >= min.events, ])
      dels.hypotheses <- temp.dels[which(!(temp.dels %in% gene.blacklist))]
      dels.mat <- dels[dels.hypotheses,]
      
      temp.muts <- rownames(somut[apply(somut, 1, sum, na.rm = TRUE) >= min.events, ])
      muts.hypotheses <- temp.muts[which(!(temp.muts %in% gene.blacklist))]
      muts.mat <- somut[muts.hypotheses,]
      
      
      # Print info about 
      if(verbose) {
        print(paste("Removing", length(gene.blacklist),"mutSig blacklist genes from hypothesis testing"))
        
        amps.removed <- length(temp.amps) - length(amps.hypotheses)
        print(paste(length(amps.hypotheses), "useable amplifications found.", amps.removed, "removed for being on the gene blacklist."))
        
        dels.removed <- length(temp.dels) - length(dels.hypotheses)
        print(paste(length(dels.hypotheses), "useable deletions found.", dels.removed, "removed for being on the gene blacklist."))
        
        muts.removed <- length(temp.muts) - length(muts.hypotheses)
        print(paste(length(muts.hypotheses), "useable mutations found.", muts.removed, "removed for being on the gene blacklist."))
      }
      
      # Save genomic hypotheses 
      if(is.na(output.folder)){
        print("No output folder selected, saving genomic hypotheses directly to object without printing.")
      } else {
        print(paste("Writing hypotheses to:", output.folder))
        dir.create(output.folder, showWarnings = FALSE)
        write.table(amps.hypotheses, file = paste0(output.folder, "/hypotheses.amps.txt"), quote = F, sep = "\t")
        write.table(dels.hypotheses, file = paste0(output.folder, "/hypotheses.dels.txt"), quote = F, sep = "\t")
        write.table(muts.hypotheses, file = paste0(output.folder, "/hypotheses.muts.txt"), quote = F, sep = "\t")
      }
      hypotheses <<- list(mut = muts.hypotheses, del = dels.hypotheses, amp = amps.hypotheses)
      
      # do aREA association
      nes.amps <- associate.events(viper, amps.mat, min.events = min.events, event.type = "Amplifications")
      nes.dels <- associate.events(viper, dels.mat, min.events = min.events, event.type = "Deletions")
      nes.muts <- associate.events(viper, muts.mat, min.events = min.events, event.type = "Mutations")
      
      nes.fusions <- NULL
      if (!is.null(fusions)) {
          fus.hypotheses <- rownames(fusions[apply(fusions, 1, sum, na.rm = TRUE) >= min.events, ])
          hypotheses <<- list(mut = muts.hypotheses, del = dels.hypotheses, amp = amps.hypotheses, fus = fus.hypotheses)
          if(!is.na(output.folder)) {
            write.table(fus.hypotheses, file = paste0(output.folder, "/hypotheses.fusions.txt"), quote = F, sep = "\t")
          }
          nes.fusions <- associate.events(viper, fusions, min.events = min.events, event.type = "Fusions")
      }
      
      # Save aREA results if desired
      if(!is.na(output.folder)){
        save(nes.amps, nes.dels, nes.muts, nes.fusions, file = paste0(output.folder, "/aREA.rda"))
      }
      # store in the object list
      nes <<- list(amp = nes.amps, del = nes.dels, mut = nes.muts, fus = nes.fusions)
}, 

  makeInteractions = function(genomic.event.types = c("amp", "del", "mut", "fus"), cindy.only = FALSE) {
      
      # NULL TFs : from the VIPER matrix, calculate p-values of the absolute mean NES score for each.  (2-tailed test, -pnorm*2).  BH-FDR < 0.05 are sig.
      # Take everything else as the background model.
      ranks[["viper"]] <<- viper.getTFScores(viper)
      sig.tfs <- viper.getSigTFS(ranks[["viper"]])
      null.TFs <- setdiff(rownames(viper), sig.tfs)
      print(paste("Found ", length(sig.tfs), " significant TFs from VIPER scores ... "))
      print(paste("Building background model from ", length(null.TFs), " nulls..."))
      
      full.type.names <- c("Amplifications", "Deletions", "Mutations", "Fusions")
      local.interactions = list()
      for (i in seq_len(4)) {
          type <- genomic.event.types[i]
          print(paste("Performing background correction for ", full.type.names[i], " NES scores..."))
          nes.thisType <- nes[[type]]
          if (is.null(nes.thisType)) {
              next
          }
          corrected.scores <- get.diggit.empiricalQvalues(viper, nes.thisType, null.TFs)
          print("Generating final interactions...")
          local.interactions[[type]] <- sig.interactors.DIGGIT(corrected.scores, nes[[type]], pathways[["cindy"]], cindy.only = cindy.only)
      }
      
      interactions <<- local.interactions
}, 

  Rank = function(use.cindy = TRUE, genomic.event.types = c("amp", "del", "mut", "fus"), use.parallel = F, cores = 1) {
      # ranks from DIGGIT scores
      integrated.z <- list()
      for (type in genomic.event.types) {
          
          if (is.null(interactions[[type]])) {
              next
          }
          # deletions/amp CNV events need to be corrected for genome location...
          if (type == "amp" || type == "del") {
              integrated.z[[type]] <- stouffer.integrate(interactions[[type]], gene.loc.mapping)
          } else {
              integrated.z[[type]] <- stouffer.integrate(interactions[[type]], NULL)
          }
      }
      # generate integrated rankings from additional sources of evidence, including CINDy and pathway databases and/or structural databases like PrePPI
      pathway.z <- list()
      
      if(use.parallel) {
        if(cores <= 1) {
          stop("Parallel processing selected but a usable number of cores has not
                been chosen. Please enter a number > 1")
        } else {
          print(paste("Parallel processing selected, using", cores, "cores"))
        }
      }
      
      for (pathway in names(pathways)) {
          
          pathway.z[[pathway]] <- list()
          # optional: use cindy scores to generate a separate ranking
          if (pathway == "cindy" && !use.cindy) {
              next
          }
          for (type in genomic.event.types) {
              if (is.null(interactions[[type]])) {
                  next
              }
              pathway.z[[pathway]][[type]] <- pathway.diggit.intersect(interactions[[type]], 
                                                                       pathways[[pathway]], pos.nes.only = TRUE, cores)
          }
      }
      
      viper.scores <- ranks[["viper"]]
      print("Integrating all data in Bayesian conditional model...")
      ranks[["integrated"]] <<- conditional.model(viper.scores, integrated.z, pathway.z)
}, 

  Cluster = function(use.parallel = F, cores = 1) {
      
      if(use.parallel) {
        if(cores <= 1) {
          stop("Parallel processing selected but multiple number of cores have not
                  been chosen. Please enter a number > 1")
        } else {
          print(paste("Parallel processing selected, using", cores, "cores"))
        }
      }
    
      # do weighted pearson correlation, using the ranks as weights
      weights <- log(ranks[["integrated"]])^2
      weights <- weights[as.character(rownames(viper))]
      w.vipermat <- weights * viper
      print("using pearson correlation with weighted vipermat")
      dist.obj <- corDist(t(w.vipermat), method = "pearson")
      print("testing clustering options, k = 2..15")
      search.results <- clusterRange(dist.obj, range = as.numeric(c(2, 15)), step = 1, cores = cores, method = "pam")
      clustering.results <<- search.results
}, 

  saturationCalculation = function(clustering.solution = NULL, cov.fraction = 0.85) {
      
      # get clustering solution to use for calculations
      if (is.null(clustering.solution)) {
        if(is.null(sample.clustering)) {
          stop("No clustering solution provided. Provide one as an argument or save one
               to the momaObj. Quitting...")
        } else {
          clustering.solution <- sample.clustering
        }
      }
      # get coverage for each subtype
      coverage.subtypes <- list()
      tmp.summaryStats <- list()
      tmp.checkpoints <- list()
      for (clus.id in unique(clustering.solution)) {
          
          print(paste0("Analyzing cluster ", clus.id))
          viper.samples <- colnames(viper[, names(clustering.solution[clustering.solution == clus.id])])
          
          # Get subtype-specific rankings: use the main rank and include only those with high mean score
          stouffer.zscores <- apply(viper[, viper.samples], 1, function(x) {
              sum(na.omit(x))/sqrt(length(na.omit(x)))
          })
          pvals <- pnorm(sort(stouffer.zscores, decreasing = T), lower.tail = F)
          sig.active.mrs <- names(pvals[p.adjust(pvals, method = "bonferroni") < 0.01])
          # ranked list of cMRs for this subtype
          subtype.specific.MR_ranks <- sort(ranks[["integrated"]][sig.active.mrs])
          
          coverage.range <- get.coverage(.self, names(subtype.specific.MR_ranks), viper.samples, topN = 100)
          coverage.subtypes[[clus.id]] <- coverage.range
          
          # 'solve' the checkpoint for each subtype
          
          # 1) generate summary stats for mean fractional coverage
          tmp.summaryStats[[clus.id]] <- merge.genomicSaturation(coverage.range, topN = 100)
          # compute best K based on fractional coverage
          sweep <- tmp.summaryStats[[clus.id]]$fraction
          names(sweep) <- tmp.summaryStats[[clus.id]]$k
          best.k <- fit.curve.percent(sweep, frac = cov.fraction)
          # pick the top cMRs based on this
          tmp.checkpoints[[clus.id]] <- names(subtype.specific.MR_ranks[seq_len(best.k)])
      }
      genomic.saturation <<- coverage.subtypes
      coverage.summaryStats <<- tmp.summaryStats
      checkpoints <<- tmp.checkpoints
  },

  makeSaturationPlots = function(clustering.solution = NULL, important.genes = NULL, max.muts = 25, max.cnvs = 10) {
    
    # print(sample.clustering)
    # get clustering solution to use for calculations
    if (is.null(clustering.solution)) {
      if (is.null(sample.clustering)) {
        stop("No clustering solution provided. Provide one as an argument or save one
               to the momaObj. Quitting...")
      } else {
        clustering.solution <- sample.clustering
      }
    }
    
    # **** incorporate checkpoint specificity here later ****
    
    tmp.identity.plots <- list()
    
    # get subtype event tables
    # print(checkpoints)
    subtype.tables <- get.subtype.event.tables(genomic.saturation, clustering.solution, checkpoints)
    
    # get summary table of unique events added in for each regulator 
    ## could be clarified/improved
    ## also potential improvement: look for inflection points of huge jumps of new unique events and highlight those regulators in particular?
    tissue.coverage.df <- merge.data.by.subtype(genomic.saturation, clustering.solution, 100)
    
    gene2band <- gene.loc.mapping$Gene.Symbol
    names(gene2band) <- gene.loc.mapping$Cytoband
    
    # Make plots each subtype
    for (k in seq_len(length(subtype.tables))) {
      # genomic events descriptive bar plot
      samples.total <- sum(clustering.solution == k)
      print(paste0("Number of samples in cluster ", k, ": ", samples.total))
      print("Getting events to plot...")
      p.identities <- plot.events(subtype.tables[[k]], important.genes, gene2band, samples.total, max.muts = 25, max.cnv = 10)
      tmp.identity.plots[["bar.plots"]][[k]] <- p.identities
    }
    
    for (k in seq_len(length(subtype.tables))) {
      # genomic coverage plot, top 100
      #subtype.df <- tissue.coverage.df[(tissue.coverage.df$subtype == k),]
      p.coverage <- genomic.plot.small(tissue.coverage.df, fraction=0.85, tissue.cluster=k)
      tmp.identity.plots[["curve.plots"]][[k]] <- p.coverage
    }
    
    
    identity.plots <<- tmp.identity.plots
    
  }



  )

)



#' @title MOMA Constructor
#' @param viper VIPER protein activity matrix with samples as columns and rows as protein IDs
#' @param mut An indicator matrix (0/1) of mutation events with samples as columns and genes as rows
#' @param fusions An indicator matrix (0/1) of fusion events with samples as columns and genes as rows
#' @param cnv A matrix of CNV scores (typically SNP6 array scores from TCGA) with samples as columns and genes as rows
#' @param pathways A named list of lists. Each named list represents interactions between proteins (keys) and their associated partners
#' @param gene.loc.mapping A data.frame of band locations and Entrez IDs
#' @param output.folder Location to store output and intermediate results 
#' @param gene.blacklist A vector of genes to exclude from mutational/CNV/fusion analysis
#' @importFrom utils data
#' @return an instance of class momaRunner
#' @export
moma_constructor <- function(viper, mut, cnv, fusions, pathways, gene.blacklist = NA_character_, 
                             output.folder = NA_character_, gene.loc.mapping = gene.map) {
    data("gene.map")
    viper <- samplename.filter(viper)
    mut <- samplename.filter(mut)
    cnv <- samplename.filter(cnv)
    
    # validate viper matrix
    if (ncol(viper) < 2) {
        print("Error: fewer than 2 samples in VIPER matrix!")
        q()
    }
    if (nrow(viper) < 2) {
        print("Error: fewer than 2 rows in VIPER matrix!")
        q()
    }
    
    # check column overlap
    nVM <- intersect(colnames(viper), colnames(mut))
    mut <- mut[, nVM]
    print(paste("Number of samples in VIPER + Mutation data:", length(nVM)))
    if (length(nVM) < 2) {
        print("Error: VIPER and mutation matrix samples don't match!")
        q()
    }
    nVM <- intersect(colnames(viper), colnames(cnv))
    cnv <- cnv[, nVM]
    print(paste("Number of samples in VIPER + CNV data:", length(nVM)))
    if (length(nVM) < 2) {
        print("Error: VIPER and CNV matrix samples don't match!")
        q()
    }
    
    if (!is.null(fusions)) {
        nVM <- intersect(colnames(viper), colnames(fusions))
        # redo type conversion to matrix: could have a single row with this datatype so we must re-type it to guard against auto conversion to 'integer'
        fusions <- as.matrix(fusions[, nVM])
        print(paste("Number of samples in VIPER + Fusion data:", length(nVM)))
        if (length(nVM) < 2) {
            print("Error: VIPER and CNV matrix samples don't match!")
            q()
        }
    }
    
    if (!is.null(gene.loc.mapping)) {
        # verify the
        if (!is(gene.loc.mapping, "data.frame")) {
            stop("Error: gene location mapping supplied is not a data.frame!")
        } else if (!("Entrez.IDs" %in% colnames(gene.loc.mapping))) {
            stop("Error: gene location mapping supplied does not have 'Entrez.IDs' attribute!")
        } else if (!("Cytoband" %in% colnames(gene.loc.mapping))) {
            stop("Error: gene location mapping supplied does not have 'Cytoband' attribute!")
        }
        
    } else {
        print("Warning: no gene - genomic location mapping provided!")
    }
    
    
    samples.common <- intersect(intersect(colnames(viper), colnames(mut)), colnames(cnv))
    print(paste("Common samples with all data analysis: ", length(samples.common)))
    
    # check TF rows in the indexes
    for (pathway in names(pathways)) {
        print(paste("Checking labels on pathway", pathway))
        I <- intersect(rownames(viper), names(pathways[[pathway]]))
        if (length(I) == 0) {
            stop("No intersection with supplied pathway! Assuming a formatting error and quitting...")
        }
        print(paste("Found labels for ", length(I), " TFs in VIPER matrix"))
    }
    
    obj <- momaRunner$new(viper = viper, mut = mut, cnv = cnv, fusions = fusions, pathways = pathways, gene.blacklist = as.character(gene.blacklist), 
        output.folder = output.folder, gene.loc.mapping = gene.loc.mapping)
    obj
}

