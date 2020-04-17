#' @title MOMA Object 
#' @description Main class encapsulating the input data and logic of the MOMA algorithm
#' @import BiocManager
#' @import stats
#' @import qvalue
#' @import parallel
#' @import reshape2
#' @import MKmisc
#' @import RColorBrewer
#' @import methods
#' @import ggplot2
#' @import magrittr
#' @import graphics
#' @import grDevices
#' @import readr
#' @import tibble
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
Moma <- setRefClass("Moma", fields = 
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
                           clustering.results = "list", # result field
                           ranks = "list", # result field
                           hypotheses = "list", # result field
                           genomic.saturation = "list", # result field
                           coverage.summaryStats = "list", # result field
                           checkpoints = "list", # result field
                           sample.clustering = "numeric" # result field 
                      ), 
                    methods = list(
                      runDIGGIT = function(fCNV = NULL, cnvthr = 0.5, min.events = 4, verbose = FALSE) {
                        "Run DIGGIT association function to get associations for driver genomic events"
                        
                        cnv.local <- NULL
                        if (is.null(fCNV)) {
                          message("No fCNV supplied, using no CNV filter!")
                          cnv.local <- cnv
                        } else {
                          message("fCNV supplied, filtering for only functional CNVs")
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
                          message("Removing ", length(gene.blacklist),
                                  " mutSig blacklist genes from hypothesis testing")
                          
                          amps.removed <- length(temp.amps) - length(amps.hypotheses)
                          message(length(amps.hypotheses), " useable amplifications found. ", 
                                  amps.removed, " removed for being on the gene blacklist.")
                          
                          dels.removed <- length(temp.dels) - length(dels.hypotheses)
                          message(length(dels.hypotheses), " useable deletions found. ", 
                                  dels.removed, " removed for being on the gene blacklist.")
                          
                          muts.removed <- length(temp.muts) - length(muts.hypotheses)
                          message(length(muts.hypotheses), " useable mutations found. ", 
                                  muts.removed, " removed for being on the gene blacklist.")
                        }
                        
                        hypotheses <<- list(mut = muts.hypotheses, del = dels.hypotheses, amp = amps.hypotheses)
                        
                        # do aREA association
                        nes.amps <- associateEvents(viper, amps.mat, min.events = min.events, 
                                                    event.type = "Amplifications",
                                                    verbose = verbose)
                        nes.dels <- associateEvents(viper, dels.mat, min.events = min.events, 
                                                    event.type = "Deletions",
                                                    verbose = verbose)
                        nes.muts <- associateEvents(viper, muts.mat, min.events = min.events, 
                                                    event.type = "Mutations",
                                                    verbose = verbose)
                        
                        nes.fusions <- NULL
                        if (!is.na(fusions[1,1])) {
                          fus.hypotheses <- rownames(fusions[apply(fusions, 1, sum, na.rm = TRUE) >= min.events, ])
                          hypotheses <<- list(mut = muts.hypotheses, del = dels.hypotheses, 
                                              amp = amps.hypotheses, fus = fus.hypotheses)
                          if(!is.na(output.folder)) {
                            write.table(fus.hypotheses, file = paste0(output.folder, "/hypotheses.fusions.txt"), 
                                        quote = FALSE, sep = "\t")
                          }
                          nes.fusions <- associateEvents(viper, fusions, 
                                                         min.events = min.events, 
                                                         event.type = "Fusions",
                                                         verbose = verbose)
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
                        message("Found ", length(sig.tfs), " significant TFs from VIPER scores ... ")
                        message("Building background model from ", length(null.TFs), " nulls...")
                        
                        full.type.names <- c("Amplifications", "Deletions", "Mutations", "Fusions")
                        local.interactions = list()
                        for (i in seq_len(4)) {
                          type <- genomic.event.types[i]
                          message("Performing background correction for ", 
                                  full.type.names[i], " NES scores...")
                          nes.thisType <- nes[[type]]
                          if (is.null(nes.thisType)) {
                            next
                          }
                          corrected.scores <- getDiggitEmpiricalQvalues(viper, nes.thisType, 
                                                                        null.TFs)
                          message("Generating final interactions...")
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
                            message("Parallel processing selected, using", cores, "cores")
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
                        message("Integrating all data in Bayesian conditional model...")
                        ranks[["integrated"]] <<- conditionalModel(viper.scores, 
                                                                   integrated.z, 
                                                                   pathway.z)
                      }, 
                      
                      Cluster = function(clus.eval = c("reliability", "silhouette"), use.parallel = FALSE, cores = 1) {
                        "Cluster the samples after applying the MOMA weights to the VIPER scores"
                        
                        if(use.parallel) {
                          if(cores <= 1) {
                            stop("Parallel processing selected but multiple number of cores have not
                  been chosen. Please enter a number > 1")
                          } else {
                            message("Parallel processing selected, using ", cores, " cores")
                          }
                        }
                        
                        # do weighted pearson correlation, using the ranks as weights
                        weights <- log(ranks[["integrated"]])^2
                        weights <- weights[as.character(rownames(viper))]
                        w.vipermat <- weights * viper
                        message("using pearson correlation with weighted vipermat")
                        dist.obj <- MKmisc::corDist(t(w.vipermat), method = "pearson")
                        message("testing clustering options, k = 2..15")
                        search.results <- clusterRange(dist.obj, range = as.numeric(c(2, 15)), 
                                                       step = 1, cores = cores, method = "pam")
                        clustering.results <<- search.results
                        
                        # save the solution with the highest reliability/silhouette to the object
                        clus.eval <- match.arg(clus.eval)
                        if(clus.eval == "reliability"){
                          top.sol <- which.max(search.results$all.cluster.reliability)
                          sample.clustering <<- search.results[[top.sol]]$clustering
                        } else if (clus.eval == "silhouette") {
                          top.sol <- which.max(search.results$all.sil.avgs)
                          sample.clustering <<- search.results[[top.sol]]$clustering
                        } else {
                          stop('Invalid clustering evaluation method provide. Choose either "reliability" or "silhouette"')
                        }
                        
                        message("Using ", clus.eval, " scores to select best solution.")
                        message("Best solution is: ", names(top.sol))
                        
                      }, 
                      
                      saturationCalculation = function(clustering.solution = NULL, cov.fraction = 0.85, 
                                                       topN = 100, verbose = FALSE) {
                        "Calculate the number of MRs it takes to represent the desired coverage fraction of events"
                        
                        # get clustering solution to use for calculations
                        if (is.null(clustering.solution)) {
                          if(is.null(sample.clustering)) {
                            stop("No clustering solution provided. Provide one as an argument or save one to the Moma Object. Quitting...")
                          } else {
                            clustering.solution <- sample.clustering
                          }
                        }
                        
                        # make sure submitted clustering solution has sufficiently large clusters (at least 5 samples)
                        clus.sizes <- table(clustering.solution) < 5
                        
                        if(any(clus.sizes)) {
                          stop("At least one cluster does not have sufficient samples. Select a different clustering solution and resave to object")
                        }
                        
                        # get coverage for each subtype
                        coverage.subtypes <- list()
                        tmp.summaryStats <- list()
                        tmp.checkpoints <- list()
                        for (clus.id in unique(clustering.solution)) {
                          
                          message("Analyzing cluster ", clus.id, " coverage using the top ", 
                                  topN, " regulators")
                          viper.samples <- colnames(viper[, names(clustering.solution[clustering.solution == clus.id])])
                          
                          # Get subtype-specific rankings: use the main rank and include only those with high mean score
                          stouffer.zscores <- apply(viper[, viper.samples], 1, function(x) {
                            sum(na.omit(x))/sqrt(length(na.omit(x)))
                          })
                          pvals <- pnorm(sort(stouffer.zscores, decreasing = TRUE), lower.tail = FALSE)
                          sig.active.mrs <- names(pvals[p.adjust(pvals, method = "bonferroni") < 0.01])
                          
                          # ranked list of cMRs for this subtype
                          subtype.specific.MR_ranks <- sort(ranks[["integrated"]][sig.active.mrs])
                          
                          # calculate per sample event coverage
                          coverage.range <- getCoverage(.self, names(subtype.specific.MR_ranks), 
                                                        viper.samples, topN = topN, verbose = verbose)
                          coverage.subtypes[[clus.id]] <- coverage.range
                          
                          # 'solve' the checkpoint for each subtype
                          
                          # 1) generate summary stats for mean fractional coverage
                          tmp.summaryStats[[clus.id]] <- mergeGenomicSaturation(coverage.range, 
                                                                                topN = topN)
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
                      },
                      saveData = function(.self, output.folder, ...) {
                        inputs <- unlist(list(...))
                        if(length(inputs) == 0){
                          message("No specific data selected to save. Saving all...")
                          to.save <- c("nes", "interactions", "clustering.results",
                                       "ranks", "hypotheses", "genomic.saturation",
                                       "coverage.summaryStats", "checkpoints")
                        } else {
                          to.save <- intersect(inputs, c("nes", "interactions", "clustering.results",
                                                         "ranks", "hypotheses", "genomic.saturation",
                                                         "coverage.summaryStats", "checkpoints"))
                          
                          if(length(to.save) == 0){
                            stop("Incorrect names supplied. Make sure names match the names of the fields in the Moma Class.")
                          } else {
                            message("Saving the following: ", paste(to.save, collapse = " "))
                          }
                          
                        }
                        
                        # check if directory exists if not create it
                        if(!dir.exists(output.folder)){
                          dir.create(output.folder)
                        }
                        
                        # create files 
                        options(max.print= 9999999)
                        options(width = 10000)
                        for(name in to.save) {
                          # if saving nes matrices save them as a table
                          if(name == "nes"){
                            for (type in names(.self[["nes"]])) {
                              write.table(.self[["nes"]][[type]], 
                                          file = paste0(output.folder, "/", type, ".nes.txt"), 
                                          quote = FALSE, sep = "\t", col.names = NA)
                            }
                          }
                          
                          # everything else save via sink
                          sink(file = paste0(output.folder, "/", name, ".txt"))
                          print(.self[[name]])
                          sink(file = NULL)
                          
                        }
                      }
                      
                    )
                    
)

utils::globalVariables(c("gene.map"))

#' MOMA Constructor Function
#' 
#' Create MOMA Object from either a MultiAssayExperiment object or a list of
#' assays.
#' See vignette for more information on how to set up and run the MOMA object
#' @param x A MultiAssayExerperiment object or list object with the following assays:
#' (note: by default assays must have these exact names. Otherwise they can be changed
#' using the viperAssay, mutMat, cnvMat and fusionMat parameters.)
#' \describe{
#' \item{viper}{VIPER protein activity matrix with samples as columns 
#' and rows as protein IDs}
#' \item{mut}{An indicator matrix (0/1) of mutation events with samples as 
#' columns and genes as rows}
#' \item{cnv}{A matrix of CNV scores (typically SNP6 array scores from TCGA) 
#' with samples as columns and genes as rows}
#' \item{fusion}{An indicator matrix (0/1) of fusion events with samples as 
#' columns and genes as rows} }
#' @param pathways A named list of lists. Each named list represents 
#' interactions between proteins (keys) and their associated partners
#' @param gene.loc.mapping A data.frame of band locations and Entrez IDs
#' @param output.folder Location to store output and intermediate results 
#' @param gene.blacklist A vector of genes to exclude from the analysis
#' @param viperAssay name associated with the viper assay in the assay object
#' @param mutMat name associated with the mutation matrix in the assay object
#' @param cnvMat name associated with the cnv matrix in the assay object
#' @param fusionMat name associated with the fusion matrix in the assay object
#' @importFrom utils data
#' @importFrom MultiAssayExperiment assays colData intersectColumns
#' @examples 
#' momaObj <- MomaConstructor(example.gbm.mae, gbm.pathways)
#' @return an instance of class Moma
#' @export
MomaConstructor <- function(x, pathways, gene.blacklist = NA_character_, 
                            output.folder = NA_character_, 
                            gene.loc.mapping = gene.map, viperAssay = "viper",
                            mutMat = "mut", cnvMat = "cnv", fusionMat = "fusion"){
  
  utils::data("gene.map")
  
  # check if x is MultiAssayExperiment or list object with the required assays
  if(is(x, "MultiAssayExperiment")){
    x <- checkMAE(x)
    type <- "mae"
  } else if (is.list(x)){
    x <- checkList(x)
    type <- "assaylist"
  } else {
    stop("Object supplied is not a MultiAssayExperiment or list object with assays.")
  }
  
  # check that genemap is exists/is valid if a different one is supplied
  checkGeneMap(gene.loc.mapping)
  
  # check pathway formatting and names
  checkPathways(pathways, x, type)
  
  ##### initialize new instance of class Moma ####
  if(type == "mae") {
    
    # first check for fusions
    if(fusionMat %in% names(assays(x))) {
      fusion <- assays(x)[[fusionMat]]
    } else {
      fusion <- matrix(NA)
    }
    
    obj <- Moma$new(viper = assays(x)[[viperAssay]], mut = assays(x)[[mutMat]], 
                    cnv = assays(x)[[cnvMat]], 
                    fusions = fusion, pathways = pathways, 
                    gene.blacklist = as.character(gene.blacklist), 
                    output.folder = output.folder, 
                    gene.loc.mapping = gene.loc.mapping)
  
    } else if (type == "assaylist") {
    
    # first check for fusions
    if(fusionMat %in% names(x)) {
      fusion <- x[[fusionMat]]
    } else {
      fusion <- matrix(NA)
    }
    
    obj <- Moma$new(viper = x[[viperAssay]], mut = x[[mutMat]], 
                    cnv = x[[cnvMat]], 
                    fusions = fusion, pathways = pathways, 
                    gene.blacklist = as.character(gene.blacklist), 
                    output.folder = output.folder, 
                    gene.loc.mapping = gene.loc.mapping)
  }
  
  obj
}



#' Check MultiAssayExperiment
#' 
#' @param mae MultiAssayExperiment object
#' @importFrom MultiAssayExperiment assays colData intersectColumns
#' @return updated/filtered MAE
#' @keywords internal
checkMAE <- function(mae){
  # confirm mae has viper, cnv, and mut assays (and optional fusion)
  if(length(intersect(names(assays(mae)) , c("viper", "mut", "cnv"))) == 3){
    message("Found the following assays:", paste(names(assays(mae)), collapse = ", "))
  } else {
    stop("Necessary assays not found. Please ensure names are 'viper', 'mut', 'cnv' (and 'fusion' if available) for the respective assays")
  }
  
  # confirm there are a sufficient number of samples
  if(nrow(colData(mae)) < 2) {
    stop("Not enough samples")
  }
  
  # print number of samples found in overlap
  message("Common samples across main assays: ", 
          suppressMessages(nrow(colData(suppressMessages(intersectColumns(mae[,,c("viper", "cnv", "mut")]))))))
  
  # filter mae to only return samples that are in the viper matrix
  vipersamples <- colnames(assays(mae)$viper)
  mae <- mae[ , vipersamples]
  mae
  
}

#' Check List of Assays
#' 
#' @param assaylist list of assays (viper, cnv, mut and fusion)
#' @return updated/filter assaylist obj
#' @keywords internal
checkList <- function(assaylist){
  if(length(intersect(names(assaylist) , c("viper", "mut", "cnv"))) == 3){
    message("Found the following assays:", paste(names(assaylist), collapse = ", "))
  } else {
    stop("Necessary assays not found. Please ensure names are 'viper', 'mut', 'cnv' (and 'fusion' if available) for the respective assays")
  }
  
  ### Confirm that a valid viper matrix has been supplied
  if (!is.matrix(assaylist$viper) || ncol(assaylist$viper) < 2 || nrow(assaylist$viper) < 2) {
    stop("Too few samples or TFs in viper matrix. Please supply valid matrix")
  }
  
  ### check for sample overlap 
  nVM <- intersect(colnames(assaylist$viper), colnames(assaylist$mut))
  assaylist$mut <- assaylist$mut[, nVM, drop = FALSE]
  message("Number of samples in VIPER + Mutation data: ", length(nVM))
  if (length(nVM) < 2) {
    stop("VIPER and mutation matrix samples don't match!")
  }
  nVM <- intersect(colnames(assaylist$viper), colnames(assaylist$cnv))
  assaylist$cnv <- assaylist$cnv[, nVM, drop = FALSE]
  message("Number of samples in VIPER + CNV data: ", length(nVM))
  if (length(nVM) < 2) {
    stop("VIPER and CNV matrix samples don't match!")
  }
  
  if ("fusion" %in% names(assaylist)) {
    nVM <- intersect(colnames(assaylist$viper), colnames(assaylist$fusion))
    # redo type conversion to matrix: could have a single row with this 
    # datatype so we must re-type it to guard against auto conversion to 'integer'
    assaylist$fusion <- as.matrix(assaylist$fusion[, nVM, drop = FALSE])
    message("Number of samples in VIPER + Fusion data: ", length(nVM))
    if (length(nVM) < 2) {
      stop("VIPER and Fusions matrix samples don't match!")
    }
  } 
  assaylist
  
}

#' Check Gene Map
#' 
#' @param gene.loc.mapping dataframe with gene names, entrez ids and cytoband locations
#' @return nothing
#' @keywords internal
checkGeneMap <- function(gene.loc.mapping){
  if (!is.null(gene.loc.mapping)) {
    # verify the column names
    if (!is(gene.loc.mapping, "data.frame")) {
      stop("Gene location mapping supplied is not a data.frame!")
    } else if (!("Entrez.IDs" %in% colnames(gene.loc.mapping))) {
      stop("Gene location mapping supplied does not have 'Entrez.IDs' attribute!")
    } else if (!("Cytoband" %in% colnames(gene.loc.mapping))) {
      stop("Gene location mapping supplied does not have 'Cytoband' attribute!")
    } else if (!("Gene.Symbol" %in% colnames(gene.loc.mapping))){
      stop("Gene location mapping supplied does not have 'Gene.Symbol' attribute!")
    }
  } else {
    stop("No gene - genomic location mapping provided!")
  }
}

#' Check Pathways
#' 
#' @param pathways A named list of lists. Each named list represents 
#' interactions between proteins (keys) and their associated partners
#' @param x the MAE or Assaylist 
#' @param type whether x is MAE or Assaylist
#' @return nothing
#' @keywords internal
checkPathways <- function(pathways, x, type){
  if(type == "mae") {
    viperTFs <- rownames(assays(x)$viper)
  } else if (type == "assaylist") {
    viperTFs <- rownames(x$viper)
  }
  
  ### Check overlap in viper protein names and pathways
  for (pathway in names(pathways)) {
    message("Checking labels on pathway ", pathway)
    I <- intersect(viperTFs, names(pathways[[pathway]]))
    if (length(I) == 0) {
      stop("No intersection with supplied pathway! Double check formatting")
    }
    message("Found labels for ", length(I), " TFs in VIPER matrix")
  }
}




#' Retain TCGA sample ids without the final letter designation ('A/B/C') 
#' 
#' @param input Matrix of expression or protein activity scores. Columns are 
#' sample names, rows are genes. Input can also just be an input vector of 
#' sample names.
#' @param desired.len length to reduce strings to. Default is 15 because of 
#' TCGA naming conventions
#' @examples 
#' sample.names <- c("TCGA-14-1825-01A", "TCGA-76-4931-01B", "TCGA-06-5418-01A")
#' sampleNameFilter(sample.names)
#' @return An identical matrix with new (shorter) column names, 
#' or a vector with the shortened names. 
#' @export
sampleNameFilter <- function(input, desired.len = 15) {
  # filter down to sample Id without the 'A/B/C sample class'.
  if(is.matrix(input)) {
    sample.ids <- vapply(colnames(input), function(x) substr(x, 1, desired.len), 
                         FUN.VALUE = character(1))
    colnames(input) <- sample.ids
  } else if (is.vector(input)) {
    input <- substr(input, 1, desired.len)
  }
  
  input
}

