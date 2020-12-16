#' @title MOMA Object 
#' @description Main class encapsulating the input data and logic of the MOMA algorithm
#' @import ggplot2
#' @import graphics
#' @import grDevices
#' @import magrittr
#' @import methods
#' @import MKmisc
#' @import parallel
#' @import poolr
#' @import qvalue
#' @import RColorBrewer
#' @import readr
#' @import reshape2
#' @import stats
#' @import tibble
#' @importFrom diggit mutualInfo correlation
#' @importFrom dplyr arrange mutate select everything bind_rows left_join across group_by summarize ungroup rename distinct
#' @importFrom parallel mc.reset.stream
#' @importFrom pracma gradient
#' @importFrom purrr map
#' @importFrom rlang .data
#' @importFrom stringr str_sub
#' @importFrom tidyr replace_na
#' @field viper matrix of inferred activity score inferred by viper
#' @field mut binary mutation matrix 1 for presence of mutation, 0 for not, NA 
#' if not determined
#' @field cnv matrix of cnv values. Can be binary or a range. 
#' @field fusions binary matrix of fusion events if applicable
#' @field expression gene expression matrix, used if fCNV or viper analyses need to be done
#' @field pathways list of pathways/connections to consider as extra evidence 
#' in the analysis
#' @field gene.blacklist character vector of genes to not include because of 
#' high mutation frequency
#' @field output.folder character vector of location to save files if desired
#' @field gene.loc.mapping data frame of gene names, entrez ids and cytoband locations
#' @field nes field for saving Normalized Enrichment Matrices from the associate events step
#' @field hypotheses results field for saving events that have enough occurences to be considered 
#' @field interactions field for saving the MR-interactions list
#' @field interactions.new field for saving updated interactions result
#' @field interactions.byCluster field for saving interactions results on a cluster specific basis
#' @field clustering.results results from clustering are saved here
#' @field ranks results field for ranking of MRs based on event association analysis
#' @field ranks.new field for new ranking results
#' @field ranks.byCluster field for saving cluster specific rank results
#' @field genomic.saturation results field for genomic saturation analysis
#' @field genomic.saturation.byCluster result field for new genomic saturation results
#' @field coverage.summaryStats results field for genomic saturation analysis
#' @field coverage.summaryStats.byCluster results field for new coverage analysis
#' @field checkpoints results field with the MRs determined to be the checkpoint for each cluster
#' @field checkpoints.byCluster resulst field for new checkpoints analysis
#' @field sample.clustering field to save sample clustering vector. Numbers are 
#' cluster assignments, names are sample ids
#' @export
Moma <- setRefClass("Moma", fields = 
                      list(viper = "matrix", 
                           mut = "matrix", 
                           cnv = "matrix", 
                           fusions = "matrix",
                           expression = "matrix",
                           pathways = "list", 
                           gene.blacklist = "character",
                           output.folder = "character", 
                           gene.loc.mapping = "data.frame",
                           fCNVs = "character", 
                           nes = "list", # result field
                           hypotheses = "list", # result field
                           interactions = "list", # result field
                           interactions.new = "data.frame", # NEW result field
                           interactions.byCluster = "list", # NEW result field
                           clustering.results = "list", # result field
                           ranks = "list", # result field
                           ranks.new = "data.frame", # NEW result field
                           ranks.byCluster = "list", # NEW result field
                           genomic.saturation = "list", # result field
                           genomic.saturation.byCluster = "list", # NEW result field
                           coverage.summaryStats = "list", # result field
                           coverage.summaryStats.byCluster = "list", # NEW result field
                           checkpoints = "list", # result field
                           checkpoints.byCluster = "list", # NEW result field
                           sample.clustering = "numeric" # result field 
                      ), 
                    methods = list(
                      runDIGGIT = function(fCNV = NULL, fCNV.sig = 0.05, fCNV.method = "mi", expression = NULL, cnvthr = 0.5, min.events = 4, verbose = FALSE) {
                        "Run DIGGIT association function to get associations for driver genomic events"
                        
                        # perform functional CNV analysis to filter out CNVs 
                        # that do not affect expression of that gene
                        
                        cnv.local <- NULL
                        if (is.null(fCNV) | isFALSE(fCNV)) {
                          message("No fCNV information supplied, using no filter to select for functional CNVs!")
                          cnv.local <- cnv
                        } else if (isTRUE(fCNV)) {
                          message("Running fCNV on supplied cnv matrix. Functional CNVs will be determined significant at a threshold of p < ", fCNV.sig)
                          
                          fCNV.tmp <- fcnvAnalysis(.self, fCNV.sig, fCNV.method, expression)
                          cnv.local <- cnv[intersect(fCNV.tmp, rownames(cnv)), ]
                          message("fCNV analysis successfully run.")
                          fCNVs <<- fCNV.tmp
                          
                        } else {
                          message("fCNV supplied, filtering for only functional CNVs")
                          cnv.local <- cnv[intersect(fCNV, rownames(cnv)), ]
                          fCNVs <<- fCNV
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
                        
                        events <- list(mut = muts.hypotheses, del = dels.hypotheses, amp = amps.hypotheses)
                        matrices <- list(mut = muts.mat, del = dels.mat, amp = amps.mat)
                        hypotheses <<- list(events = events, matrices = matrices)
                        
                        # hypotheses <<- list(mut = muts.hypotheses, del = dels.hypotheses, amp = amps.hypotheses,
                        #                     mut.mat = muts.mat, del.mat = dels.mat, amp.mat = amps.mat)
                        
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
                          fus.mat <- fusions[fus.hypotheses,]
                          
                          events <- list(mut = muts.hypotheses, del = dels.hypotheses, amp = amps.hypotheses, fus = fus.hypotheses)
                          matrices <- list(mut = muts.mat, del = dels.mat, amp = amps.mat, fus = fus.mat)
                          hypotheses <<- list(events = events, matrices = matrices)
                          
                          # hypotheses <<- list(mut = muts.hypotheses, del = dels.hypotheses, 
                          #                     amp = amps.hypotheses, fus = fus.hypotheses)
                          nes.fusions <- associateEvents(viper, fusions, 
                                                         min.events = min.events, 
                                                         event.type = "Fusions",
                                                         verbose = verbose)
                        }
                        
                        # store results in the object list
                        nes <<- list(amp = nes.amps, del = nes.dels, mut = nes.muts, fus = nes.fusions)
                        
                        
                        ### TODO: Add in new DIGGIT null model here potentially?
                        
                        
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
                      
                      makeInteractionsNew = function(genomic.event.types = c("amp", "del", "mut", "fus")){
                        "Make merged file of each MR-Event pairing with all associated pathway test values"
                        
                        # initiate tibble for results
                        full.interaction.table <- tibble::tibble(.rows = 0)
                        
                        # for each event type gather nes/aqtl scores, and whatever other pathway scores are available
                        for(type in genomic.event.types){
                          
                          # matrix of nes/aqtl scores. columns are TFs and rows are events
                          nes.thisType <- nes[[type]]
                          
                          if(is.null(nes.thisType)) next
                          
                          message("Getting interactions for ", type, " events...")
                          
                          interactions.df <- reshape2::melt(nes.thisType, varnames = c("event", "regulator"), value.name = "aQTL") %>%
                            dplyr::select(regulator, dplyr::everything()) %>% 
                            dplyr::mutate(dplyr::across(.cols = c(regulator, event), as.character))
                          
                          # go through pathway lists to find matching data (if it exists)
                          # TODO: add in option where this section checks if the data is list or dataframe
                          # in order to skip the unlisting steps
                          for(p in seq_along(pathways)) {
                            
                            p.name <- names(pathways)[[p]]
                            pathway <- pathways[[p]]
                            
                            # convert named list of MR-Events to a dataframe
                            full.df <- tibble::tibble(.rows = 0)
                            for(mr in names(pathway)) {
                              df <- tibble::enframe(pathway[[mr]], 
                                                    name = "event", value = "pval") %>%
                                dplyr::mutate(regulator = mr)
                              full.df <- dplyr::bind_rows(full.df, df)
                              
                            }
                            
                            # merge the pathway dataframe to interaction.df
                            # change pval to pathway name so it doesn't get written over
                            interactions.df <- dplyr::left_join(interactions.df, full.df, by = c("regulator", "event")) %>%
                              dplyr::rename(!!p.name := pval)
                            
                          }
                          
                          # add event type to table then merge with full.interaction.table
                          interactions.df$type <- type
                          full.interaction.table <- dplyr::bind_rows(full.interaction.table, interactions.df) 
                          
                        }
                        
                        # only keep the rows that are distinct 
                        # (preppi interaction data can incur duplicates)
                        # save aQTL direction information as a new column then convert to p-values
                        full.interaction.table <- full.interaction.table %>% 
                          dplyr::distinct() %>%
                          dplyr::mutate(sign = sign(aQTL),
                                        aQTL = 2*pnorm(-abs(aQTL)))
                        
                       
                        interactions.new <<- full.interaction.table
                        
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
                      
                      RankNew = function(na.value = NA, aQTL.threshold = 1) {
                        "Combine all genomic information to create event-MR ranking and full MR ranking"
                        
                        ## convert all scores to newly normalized p-values
                        # filter to events that are below aQTL threshold
                        interactions.temp <- rankNormalize(interactions.new) %>%
                          dplyr::filter(aQTL <= aQTL.threshold)
                        
                        
                        ## Two part integration of pvalues 
                        # first across MR-Event pairs
                        # then all events associated with an MR
                        message("Rank integrating regulator events...")
                        res <- twoStepRankIntegration(interactions.temp, na.value)
                        
                        # save to main object
                        interactions.new <<- res$interactions.df
                        ranks.new <<- res$ranks.df
                        
                      },
                      
                      Cluster = function(clus.eval = c("reliability", "silhouette"), use.parallel = FALSE, cores = 1, new = FALSE) {
                        "Cluster the samples after applying the MOMA weights to the VIPER scores"
                        
                        ## TODO: Integrate iterative clustering? see about updating iterClust directly then including it
                        ## TODO: Try different ways of handling weighting 
                        ##       - log(pvals)? adjusted pvals? just ranks?
                        
                        if(use.parallel) {
                          if(cores <= 1) {
                            stop("Parallel processing selected but multiple number of cores have not
                  been chosen. Please enter a number > 1")
                          } else {
                            message("Parallel processing selected, using ", cores, " cores")
                          }
                        }
                        
                        # do weighted pearson correlation, using the ranks as weights
                        
                        # make option for new/old run (for vignette issue)
                        if(isFALSE(new)) {
                          weights <- log(ranks[["integrated"]])^2
                        } else {
                          weights <- tibble::deframe(ranks.new)
                          weights <- weights[as.character(rownames(viper))]
                          
                          # adjustment of weights
                          # TODO 
                          # currently just taking log. subject to change
                          
                          weights <- -log(weights)
                        }
                        
                        
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
                      
                      RankByCluster = function(clustering.solution = NULL, min.events.per.cluster = 2, 
                                               na.value = NA, aQTL.threshold = 1) {
                        "Recalculate regulator rankings for each cluster based on presence of genomic events in that cluster"
                        
                        # get clustering solution to use for calculations
                        if (is.null(clustering.solution)) {
                          if (is.null(sample.clustering)) {
                            stop("No clustering solution provided. Provide one as an argument or save one
                                 to the momaObj. Quitting...")
                          } else {
                            clustering.solution <- sample.clustering
                          }
                        }
                        
                        # get previously filtered gene hypotheses from runDIGGIT
                        if(is.null(hypotheses)) {
                          stop("No hypotheses found. Do runDIGGIT() before continuing.")
                        }
                        
                        event.types <- unique(interactions.new$type)
                        
                        ranks.byCluster.tmp <- list()
                        interactions.byCluster.tmp <- list()
                        for (cluster in unique(clustering.solution)) {
                          message(paste0("Integrating ranks for cluster ", cluster))
                          
                          inCluster.samples <- names(clustering.solution[clustering.solution==cluster])
                          
                          # genomic events specific to this cluster: do over representation analysis
                          # filter interactions list to only these events
                          overrep.events <- list()
                          filtered.interactions <- tibble::tibble(.rows = 0)
                          for(etype in event.types) {
                            res <- overrep.analysis(hypotheses[["matrices"]][[etype]], inCluster.samples, min.events.per.cluster)
                            overrep.events[[etype]] <- res
                            
                            res.interactions <- interactions.new %>% 
                              dplyr::filter(type == etype) %>%
                              dplyr::filter(event %in% res)
                            
                            filtered.interactions <- dplyr::bind_rows(filtered.interactions, res.interactions)
                          }
                          
                          filtered.interactions$cluster <- cluster
                          
                          # integrate events to get a rank for each regulator
                          ranks.df <- filtered.interactions %>% dplyr::group_by(regulator) %>%
                            dplyr::summarize(int.mr.p = poolr::fisher(int.p)$p) %>% dplyr::arrange(int.mr.p)
                          
                          ranks.df$int.mr.p <- stats::p.adjust(ranks.df$int.mr.p, method = "fdr")
                          ranks.df$cluster <- cluster
                          
                          ranks.byCluster.tmp[[cluster]] <- ranks.df
                          interactions.byCluster.tmp[[cluster]] <- filtered.interactions
                          
                        }
                        
                        # save to main object
                        ranks.byCluster <<- ranks.byCluster.tmp
                        interactions.byCluster <<- interactions.byCluster.tmp
                        
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
                        checkpoints <<- setNames(tmp.checkpoints, seq_along(tmp.checkpoints))
                      },
                      saturationCalculationNew = function(clustering.solution = NULL, cytoband.collapse = T, topN = 100) {
                        "Calculate the number of MRs it takes to represent the desired coverage fraction of events"
                        
                        # main update from old version: 
                        # calculation of saturation point as when the derivative of the curve is 0
                        # ie the addition of more regulators does not add events
                        
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
                          message("Analyzing cluster ", clus.id, " coverage...")
                          
                          viper.samples <- colnames(viper[, names(clustering.solution[clustering.solution == clus.id])])
                          
                          # Get subtype-specific rankings: use the main rank and include only those with significantly high/low mean score
                          stouffer.zscores <- apply(viper[, viper.samples], 1, function(x) {
                            sum(na.omit(x))/sqrt(length(na.omit(x)))
                          })
                          
                          pvals <- sort(2*pnorm(-abs(stouffer.zscores)), decreasing = F)
                          sig.active.mrs <- names(pvals[p.adjust(pvals, method='bonferroni') < 0.25])
                          
                          # rank using the subtype-specific rankings generated in this function, above. 
                          # otherwise this is the same analysis done on the overall rankings
                          subtype.specific.MR_ranks <- ranks.byCluster[[clus.id]] %>% 
                            dplyr::filter(regulator %in% sig.active.mrs) %>%
                            dplyr::select(regulator) %>% unlist(use.names = F)
                          
                          
                          # filter cluster interactions list to only the significant MRs and significant events
                          clus.interactions.topmrs <- interactions.byCluster[[clus.id]] %>%
                            dplyr::filter(regulator %in% subtype.specific.MR_ranks & int.p < 0.05) 
                          
                          clus.interactions.topmrs <- dplyr::left_join(clus.interactions.topmrs, gene.loc.mapping, 
                                                                       by = c("event" = "Entrez.IDs"))
                          
                          # if(isTRUE(cytoband.collapse)) {
                          #   # merge in gene map so that close amps/dels can be consolidated to one event
                          #   
                          #   
                          #   cnvs.cytoband.collapse <- clus.interactions.topmrs %>% 
                          #     dplyr::filter(type %in% c("amp", "del")) %>%
                          #     dplyr::group_by(regulator, type, sign, Cytoband) %>% 
                          #     dplyr::summarize(int.p = min(int.p)) %>% 
                          #     dplyr::ungroup() %>% 
                          #     dplyr::left_join(clus.interactions.topmrs)
                          #   
                          #   other.interactions <- clus.interactions.topmrs %>%
                          #     dplyr::filter(!type %in% c("amp", "del"))
                          #   
                          #   clus.interactions.topmrs <- dplyr::bind_rows(cnvs.cytoband.collapse, other.interactions)
                          #   
                          # }
                          
                          coverage.range <- sampleEventOverlap(.self, viper.samples, subtype.specific.MR_ranks, clus.interactions.topmrs,
                                                               cytoband.collapse, topN)
                          coverage.subtypes[[clus.id]] <- coverage.range
                          
                          # Solve the checkpoint for each subtype
                          # Merge information about each sample together
                          # Then use discrete numerical gradient function to calculate the derivative at each point
                          # of the fraction function. Choose first instance of 0 / inflection
                          # Pick top cMRs based on this
                          tmp.summaryStats[[clus.id]] <- genomicSaturationSummary(coverage.range, topN)
                          
                          best.k <- getInflection(tmp.summaryStats[[clus.id]]$fraction, clus.id)
                          
                          tmp.checkpoints[[clus.id]] <- subtype.specific.MR_ranks[seq_len(best.k)]
                          }
                        
                        checkpoints.byCluster <<- tmp.checkpoints
                        genomic.saturation.byCluster <<- coverage.subtypes
                        coverage.summaryStats.byCluster <<- tmp.summaryStats
                        
                      },
                      saveData = function(.self, output.folder, ...) {
                        inputs <- unlist(list(...))
                        if(length(inputs) == 0){
                          message("No specific data selected to save. Saving all...")
                          to.save <- c("nes", "interactions", "clustering.results",
                                       "sample.clustering","ranks", "hypotheses", 
                                       "genomic.saturation","coverage.summaryStats", "checkpoints")
                        } else {
                          to.save <- intersect(inputs, c("nes", "interactions", "clustering.results", "sample.clustering",
                                                         "ranks", "hypotheses", "genomic.saturation",
                                                         "coverage.summaryStats", "checkpoints"))
                          
                          if(length(to.save) == 0){
                            stop("Incorrect names supplied. Make sure names match the names of the results fields in the Moma Class.")
                          } else {
                            message("Saving the following: ", paste(to.save, collapse = " "))
                          }
                          
                        }
                        
                        # check if output.folder name ends in / or \\ or not
                        if(!stringr::str_sub(output.folder, start = -1) %in% c("/", "\\")){
                          stop("Output folder does not end with a slash")
                        }
                        
                        # check if directory exists if not create it
                        if(!dir.exists(output.folder)){
                          dir.create(output.folder)
                        }
                        
                        ### create files
                        for(name in to.save) {
                          
                          # confirm that result exists before trying to save
                          if(length(.self[[name]]) == 0) {
                            message("No results found for ", name, ". Skipping...")
                            next
                          }
                          
                          # save nes matrices as tables
                          if(name == "nes"){
                            for (type in names(.self[[name]])) {
                              write.table(.self[[name]][[type]], 
                                          file = paste0(output.folder, type, ".", name, ".txt"), 
                                          quote = FALSE, sep = "\t", col.names = NA)
                            }
                          }
                          
                          # save interactions as 3 separate dfs 
                          if(name == "interactions"){
                            for (type in names(.self[["interactions"]])){
                              full.df <- tibble::tibble(regulator = NA, event = NA, nes = NA, .rows = 0)
                              for(mr in names(.self[["interactions"]][[type]])) {
                                df <- tibble::enframe(.self[["interactions"]][[type]][[mr]], 
                                              name = "event", value = "nes") %>%
                                  dplyr::mutate(regulator = mr)
                                full.df <- dplyr::bind_rows(full.df, df)
                              }
                              
                              write.table(full.df, file = paste0(output.folder, type,".", name, ".txt"),
                                          quote = FALSE, sep = "\t", row.names = FALSE)
                            }
                          }
                          
                          # save ranks as df
                          if(name == "ranks"){
                            df <- lapply(.self[["ranks"]], tibble::enframe) %>%
                              dplyr::bind_rows() %>% 
                              dplyr::mutate(type = rep(c("viper", "integrated"), each = length(.self[["ranks"]][[1]])))
                            
                            write.table(df, file = paste0(output.folder, name, ".txt"),
                                        quote = FALSE, sep = "\t", row.names = FALSE)
                          }
                          
                          # save hypotheses as df
                          if(name == "hypotheses"){
                            df <- lapply(.self[[name]][["events"]], function(x){
                              tibble::tibble(type = NA, gene = x)
                            }) %>% dplyr::bind_rows()
                            event.nums <- vapply(.self[[name]][["events"]], length, numeric(1))
                            df$type <- rep(names(event.nums), event.nums)
                            
                            write.table(df, file = paste0(output.folder, name, ".txt"),
                                        quote = FALSE, sep = "\t", row.names = FALSE)
                            
                          }
                          
                          # save checkpoints as df
                          if(name == "checkpoints"){
                            df <- lapply(.self[[name]], function(x){
                              tibble::tibble(cluster = NA, regulator = x)
                            }) %>% dplyr::bind_rows()
                            event.nums <- vapply(.self[[name]], length, numeric(1))
                            df$cluster <- rep(names(event.nums), event.nums)
                            
                            write.table(df, file = paste0(output.folder, name, ".txt"),
                                        quote = FALSE, sep = "\t", row.names = FALSE)
                            
                          }
                          
                          # save sample clustering as df
                          if(name == "sample.clustering"){
                            df <- tibble::enframe(.self[["sample.clustering"]], 
                                                  name = "sample", value = "cluster") 
                            
                            write.table(df, file = paste0(output.folder, name, ".txt"),
                                        quote = FALSE, sep = "\t", row.names = FALSE)
                          }
                          
                          # save coverage.summary stats as df
                          if(name == "coverage.summaryStats") {
                            df <- dplyr::bind_rows(.self[[name]]) %>%
                              dplyr::mutate(cluster = rep(seq_along(.self[[name]]), each = 100))
                            
                            write.table(df, file = paste0(output.folder, name, ".txt"),
                                        quote = FALSE, sep = "\t", row.names = FALSE)
                          }
                          
                          # save clustering results and genomic saturation as lists
                          if(name %in% c("clustering.results", "genomic.saturation")){
                            results <- .self[[name]]
                            save(results, 
                                 file = paste0(output.folder, name, ".rda"))
                          }
                          
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
#' using the viperAssay, mutMat, cnvMat, fusionMat, and expressionMat parameters.)
#' \describe{
#' \item{viper}{VIPER protein activity matrix with samples as columns 
#' and rows as protein IDs}
#' \item{mut}{An indicator matrix (0/1) of mutation events with samples as 
#' columns and genes as rows}
#' \item{cnv}{A matrix of CNV scores (typically SNP6 array scores from TCGA) 
#' with samples as columns and genes as rows}
#' \item{fusion}{An indicator matrix (0/1) of fusion events with samples as 
#' columns and genes as rows} 
#' \item{expression}{A matrix of gene expression values} }
#' @param pathways A named list of lists. Each named list represents 
#' interactions between proteins (keys) and their associated partners
#' @param gene.loc.mapping A data.frame of band locations and Entrez IDs
#' @param output.folder Location to store output and intermediate results 
#' @param gene.blacklist A vector of genes to exclude from the analysis
#' @param viperAssay name associated with the viper assay in the assay object
#' @param mutMat name associated with the mutation matrix in the assay object
#' @param cnvMat name associated with the cnv matrix in the assay object
#' @param fusionMat name associated with the fusion matrix in the assay object
#' @param expressionMat name associated with the expression matrix in the assay object
#' @importFrom utils data
#' @importFrom MultiAssayExperiment assays colData intersectColumns
#' @examples 
#' momaObj <- MomaConstructor(example.gbm.mae, gbm.pathways)
#' @return an instance of class Moma
#' @export
MomaConstructor <- function(x, pathways, gene.blacklist = NA_character_, 
                            output.folder = NA_character_, 
                            gene.loc.mapping = gene.map, viperAssay = "viper",
                            mutMat = "mut", cnvMat = "cnv", fusionMat = "fusion",
                            expressionMat = "expression"){
  
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
    
    # first check for fusions and expression
    if(fusionMat %in% names(assays(x))) {
      fusion <- assays(x)[[fusionMat]]
    } else {
      fusion <- matrix(NA)
    }
    
    if(expressionMat %in% names(assays(x))) {
      exp <- assays(x)[[expressionMat]]
    } else {
      exp <- matrix(NA)
    }
    
    obj <- Moma$new(viper = assays(x)[[viperAssay]], mut = assays(x)[[mutMat]], 
                    cnv = assays(x)[[cnvMat]], 
                    fusions = fusion, expression = exp,
                    pathways = pathways, 
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
    
    if(expressionMat %in% names(x)) {
      exp <- x[[expressionMat]]
    } else {
      exp <- matrix(NA)
    }
    
    
    
    obj <- Moma$new(viper = x[[viperAssay]], mut = x[[mutMat]], 
                    cnv = x[[cnvMat]], 
                    fusions = fusion, expression = exp,
                    pathways = pathways, 
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
  
  # remove any rows from viper matrix that are all NA 
  # (happens in vipermats from full TCGA viper analysis)
  if(sum(is.na(assays(mae)$viper)) > 0) {
    na.mrs <- which(is.na(assays(mae)$viper[,1]))
    mae@ExperimentList@listData$viper <- assays(mae)$viper[-na.mrs,]
    }
  
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
  # also remove any rows that are all NAs
  if (!is.matrix(assaylist$viper) || ncol(assaylist$viper) < 2 || nrow(assaylist$viper) < 2) {
    stop("Too few samples or TFs in viper matrix. Please supply valid matrix")
  }
  
  if(sum(is.na(assaylist$viper)) > 0) {
    na.mrs <- which(is.na(assaylist$viper[,1]))
    assaylist$viper <- assaylist$viper[-na.mrs,]
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

