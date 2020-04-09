
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
#' @import moma.gbmexample
#' @field viper matrix of inferred activity score inferred by viper
#' @field mut binary mutation matrix 1 for presence of mutation, 0 for not, NA 
#' if not determined
#' @field cnv matrix of cnv values. Can be binary or a range. 
#' @field fusions binary matrix of fusion events if appliable
#' @field pathways list of pathways/connections to consider as extra evidence 
#' in the analysis
#' @field gene.blacklist character vector of genes to not include because of 
#' high mutation frequency
#' @field output.folder character vector of location to save files if desired
#' @field gene.loc.mapping data frame of gene names, entrez ids and cytoband locations
#' @field nes field for saving Normalized Enrichment Matrices from the associate events step
#' @field interactions field for saving the MR-interactions list
#' @field clustering.results results from clustering are saved here
#' @field ranks results field for ranking of MRs based on event association analysis
#' @field hypotheses results field for saving events that have enough occurences to be considered 
#' @field genomic.saturation results field for genomic saturation analysis
#' @field coverage.summaryStats results field for genomic saturation analysis
#' @field checkpoints results field with the MRs determined to be the checkpoint for each cluster
#' @field sample.clustering field to save sample clustering vector. Numbers are 
#' cluster assignments, names are sample ids
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
                                 sample.clustering = "numeric" 
                                 ), 
                          methods = list(
  runDIGGIT = function(fCNV = NULL, cnvthr = 0.5, min.events = 4, verbose = FALSE) {
    "Run DIGGIT association function to get associations for driver genomic events"
    
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
        print(paste("Removing", length(gene.blacklist),"mutSig blacklist genes 
                    from hypothesis testing"))
        
        amps.removed <- length(temp.amps) - length(amps.hypotheses)
        print(paste(length(amps.hypotheses), "useable amplifications found.", 
                    amps.removed, "removed for being on the gene blacklist."))
        
        dels.removed <- length(temp.dels) - length(dels.hypotheses)
        print(paste(length(dels.hypotheses), "useable deletions found.", 
                    dels.removed, "removed for being on the gene blacklist."))
        
        muts.removed <- length(temp.muts) - length(muts.hypotheses)
        print(paste(length(muts.hypotheses), "useable mutations found.", 
                    muts.removed, "removed for being on the gene blacklist."))
      }
      
      # Save genomic hypotheses 
      if(is.na(output.folder)){
        print("No output folder selected, saving genomic hypotheses 
              directly to object without printing.")
      } else {
        print(paste("Writing hypotheses to:", output.folder))
        dir.create(output.folder, showWarnings = FALSE)
        write.table(amps.hypotheses, file = paste0(output.folder, "/hypotheses.amps.txt"), 
                    quote = FALSE, sep = "\t")
        write.table(dels.hypotheses, file = paste0(output.folder, "/hypotheses.dels.txt"), 
                    quote = FALSE, sep = "\t")
        write.table(muts.hypotheses, file = paste0(output.folder, "/hypotheses.muts.txt"), 
                    quote = FALSE, sep = "\t")
      }
      hypotheses <<- list(mut = muts.hypotheses, del = dels.hypotheses, amp = amps.hypotheses)
      
      # do aREA association
      nes.amps <- associateEvents(viper, amps.mat, min.events = min.events, 
                                  event.type = "Amplifications")
      nes.dels <- associateEvents(viper, dels.mat, min.events = min.events, 
                                  event.type = "Deletions")
      nes.muts <- associateEvents(viper, muts.mat, min.events = min.events, 
                                  event.type = "Mutations")
      
      nes.fusions <- NULL
      if (!is.null(fusions)) {
          fus.hypotheses <- rownames(fusions[apply(fusions, 1, sum, na.rm = TRUE) >= min.events, ])
          hypotheses <<- list(mut = muts.hypotheses, del = dels.hypotheses, 
                              amp = amps.hypotheses, fus = fus.hypotheses)
          if(!is.na(output.folder)) {
            write.table(fus.hypotheses, file = paste0(output.folder, "/hypotheses.fusions.txt"), 
                        quote = FALSE, sep = "\t")
          }
          nes.fusions <- associateEvents(viper, fusions, min.events = min.events, event.type = "Fusions")
      }
      
      # Save aREA results if desired
      if(!is.na(output.folder)){
        save(nes.amps, nes.dels, nes.muts, nes.fusions, file = paste0(output.folder, "/aREA.rda"))
      }
      # store in the object list
      nes <<- list(amp = nes.amps, del = nes.dels, mut = nes.muts, fus = nes.fusions)
}, 

  makeInteractions = function(genomic.event.types = c("amp", "del", "mut", "fus"),
                              cindy.only = FALSE) {
    "Make interaction web for significant MRs based on their associated events"
    
      # NULL TFs : from the VIPER matrix, calculate p-values of the absolute mean
      # NES score for each.  (2-tailed test, -pnorm*2).  BH-FDR < 0.05 are sig.
      # Take everything else as the background model.
      ranks[["viper"]] <<- viperGetTFScores(viper)
      sig.tfs <- viperGetSigTFS(ranks[["viper"]])
      null.TFs <- setdiff(rownames(viper), sig.tfs)
      print(paste("Found ", length(sig.tfs), " significant TFs from VIPER scores ... "))
      print(paste("Building background model from ", length(null.TFs), " nulls..."))
      
      full.type.names <- c("Amplifications", "Deletions", "Mutations", "Fusions")
      local.interactions = list()
      for (i in seq_len(4)) {
          type <- genomic.event.types[i]
          print(paste("Performing background correction for ", 
                      full.type.names[i], " NES scores..."))
          nes.thisType <- nes[[type]]
          if (is.null(nes.thisType)) {
              next
          }
          corrected.scores <- getDiggitEmpiricalQvalues(viper, nes.thisType, 
                                                        null.TFs)
          print("Generating final interactions...")
          local.interactions[[type]] <- sigInteractorsDIGGIT(corrected.scores, 
                                                             nes[[type]], pathways[["cindy"]], 
                                                             cindy.only = cindy.only)
      }
      
      interactions <<- local.interactions
}, 

  Rank = function(use.cindy = TRUE, genomic.event.types = c("amp", "del", "mut", "fus"), 
                  use.parallel = FALSE, cores = 1) {
    "Rank MRs based on DIGGIT scores and number of associated events"  
    
    # ranks from DIGGIT scores
      integrated.z <- list()
      for (type in genomic.event.types) {
          
          if (is.null(interactions[[type]])) {
              next
          }
          # deletions/amp CNV events need to be corrected for genome location...
          if (type == "amp" || type == "del") {
              integrated.z[[type]] <- stoufferIntegrate(interactions[[type]], 
                                                        gene.loc.mapping)
          } else {
              integrated.z[[type]] <- stoufferIntegrate(interactions[[type]], 
                                                        NULL)
          }
      }
      # generate integrated rankings from additional sources of evidence, 
      # using CINDy and pathway databases and/or structural databases like PrePPI
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
              pathway.z[[pathway]][[type]] <- pathwayDiggitIntersect(interactions[[type]], 
                                  pathways[[pathway]], pos.nes.only = TRUE, cores)
          }
      }
      
      viper.scores <- ranks[["viper"]]
      print("Integrating all data in Bayesian conditional model...")
      ranks[["integrated"]] <<- conditionalModel(viper.scores, 
                                                 integrated.z, 
                                                 pathway.z)
}, 

  Cluster = function(use.parallel = FALSE, cores = 1) {
    "Cluster the samples after applying the MOMA weights to the VIPER scores"
      
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
      search.results <- clusterRange(dist.obj, range = as.numeric(c(2, 15)), 
                                     step = 1, cores = cores, method = "pam")
      clustering.results <<- search.results
}, 

  saturationCalculation = function(clustering.solution = NULL, cov.fraction = 0.85) {
    "Calculate the number of MRs it takes to represent the desired coverage fraction of events"
      
      # get clustering solution to use for calculations
      if (is.null(clustering.solution)) {
        if(is.null(sample.clustering)) {
          stop("No clustering solution provided. Provide one as an argument or 
                save one to the momaObj. Quitting...")
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
          pvals <- pnorm(sort(stouffer.zscores, decreasing = TRUE), lower.tail = FALSE)
          sig.active.mrs <- names(pvals[p.adjust(pvals, method = "bonferroni") < 0.01])
          # ranked list of cMRs for this subtype
          subtype.specific.MR_ranks <- sort(ranks[["integrated"]][sig.active.mrs])
          
          coverage.range <- getCoverage(.self, names(subtype.specific.MR_ranks), 
                                        viper.samples, topN = 100)
          coverage.subtypes[[clus.id]] <- coverage.range
          
          # 'solve' the checkpoint for each subtype
          
          # 1) generate summary stats for mean fractional coverage
          tmp.summaryStats[[clus.id]] <- mergeGenomicSaturation(coverage.range, 
                                                                topN = 100)
          # compute best K based on fractional coverage
          sweep <- tmp.summaryStats[[clus.id]]$fraction
          names(sweep) <- tmp.summaryStats[[clus.id]]$k
          best.k <- fitCurvePercent(sweep, frac = cov.fraction)
          # pick the top cMRs based on this
          tmp.checkpoints[[clus.id]] <- names(subtype.specific.MR_ranks[seq_len(best.k)])
      }
      genomic.saturation <<- coverage.subtypes
      coverage.summaryStats <<- tmp.summaryStats
      checkpoints <<- tmp.checkpoints
  }

  )

)

utils::globalVariables(c("gene.map"))

#' @title MOMA Constructor
#' see vignette for more information on how to set up and run the MOMA object
#' @param viper VIPER protein activity matrix with samples as columns 
#' and rows as protein IDs
#' @param mut An indicator matrix (0/1) of mutation events with samples as 
#' columns and genes as rows
#' @param fusions An indicator matrix (0/1) of fusion events with samples as 
#' columns and genes as rows
#' @param cnv A matrix of CNV scores (typically SNP6 array scores from TCGA) 
#' with samples as columns and genes as rows
#' @param pathways A named list of lists. Each named list represents 
#' interactions between proteins (keys) and their associated partners
#' @param gene.loc.mapping A data.frame of band locations and Entrez IDs
#' @param output.folder Location to store output and intermediate results 
#' @param gene.blacklist A vector of genes to exclude from the analysis
#' @importFrom utils data
#' @examples
#' library(moma.gbmexample) 
#' data("gbm.example")
#' momaObj <- momaConstructor(gbm.example$vipermat, 
#' gbm.example$rawsnp, 
#' gbm.example$rawcnv, 
#' gbm.example$fusions, 
#' pathways = list(cindy=gbm.example$cindy, preppi=gbm.example$preppi), 
#' gene.blacklist=gbm.example$mutSig)
#' 
#' @return an instance of class momaRunner
#' @export
momaConstructor <- function(viper, mut, cnv, fusions, pathways, 
                            gene.blacklist = NA_character_, 
                            output.folder = NA_character_, 
                            gene.loc.mapping = gene.map) {
    utils::data("gene.map")
    viper <- sampleNameFilter(viper)
    mut <- sampleNameFilter(mut)
    cnv <- sampleNameFilter(cnv)
    
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
        # redo type conversion to matrix: could have a single row with this 
        # datatype so we must re-type it to guard against auto conversion to 'integer'
        fusions <- as.matrix(fusions[, nVM])
        print(paste("Number of samples in VIPER + Fusion data:", length(nVM)))
        if (length(nVM) < 2) {
            print("Error: VIPER and Fusions matrix samples don't match!")
            q()
        }
    }
    
    if (!is.null(gene.loc.mapping)) {
        # verify the column names
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
            stop("No intersection with supplied pathway! 
                 Assuming a formatting error and quitting...")
        }
        print(paste("Found labels for ", length(I), " TFs in VIPER matrix"))
    }
    
    obj <- momaRunner$new(viper = viper, mut = mut, cnv = cnv, 
                          fusions = fusions, pathways = pathways, 
                          gene.blacklist = as.character(gene.blacklist), 
                          output.folder = output.folder, 
                          gene.loc.mapping = gene.loc.mapping)
    obj
}


#' Retain TCGA sample ids without the final letter designation ('A/B/C') 
#' 
#' @param input Matrix of expression or protein activity scores. Columns are 
#' sample names, rows are genes. Input can also just be an input vector of 
#' sample names.
#' @examples 
#' sample.names <- c("TCGA-14-1825-01A", "TCGA-76-4931-01B", "TCGA-06-5418-01A")
#' sampleNameFilter(sample.names)
#' @return An identical matrix with new (shorter) column names, 
#' or a vector with the shortened names. 
#' @export
sampleNameFilter <- function(input) {
  # filter down to sample Id without the 'A/B/C sample class'.
  if(is.matrix(input)) {
    sample.ids <- vapply(colnames(input), function(x) substr(x, 1, 15), 
                         FUN.VALUE = character(1))
    colnames(input) <- sample.ids
  } else if (is.vector(input)) {
    input <- substr(input, 1, 15)
  }
  
  input
}

