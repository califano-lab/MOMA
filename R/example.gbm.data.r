#' Glioblastoma (GBM) Example Dataset
#' 
#' MultiAssayExperiment Object containing all the genomic assays needed to run
#' the example code for MOMA
#' 
#' @format An MultiAssayExperiment object with 4 different sets of GBM assays
#' \describe{
#' \item{viper}{matrix of viper scores with samples in columns and regulators
#'  across the rows}
#' \item{mut}{matrix of samples and genes with potential mutations.
#' 0 for no mutation, 1 for presence of some non-silent mutation}
#' \item{cnv}{matrix of samples and genes with copy number variant scores}
#' \item{fusion}{matrix of samples and fusion genes.
#' 0 for no fusion, 1 for presence of fusion}
#' }
"example.gbm.mae"


#' Glioblastoma (GBM) Pathways
#' 
#' Object containing information about the biological pathways that will be used
#' in the analysis
#' 
#' @format A list of lists named "cindy" and "preppi" respectively
#' \describe{
#' \item{cindy}{list of regulators, each with a set of modulators and p values
#' representing their CINDY inferred association}
#' \item{preppi}{list of regulators, each with a set of potential binding
#' partners and PREPPi inferred p values for probability of binding}
#' }
"gbm.pathways"