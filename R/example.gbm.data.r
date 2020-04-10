#' Glioblastoma (GBM) Example Dataset
#'
#' Object containing all the data needed to run the example code for MOMA
#'
#' @format An object with 9 different sets of GBM data needed to run the MOMA analysis
#' \describe{
#' \item{vipermat}{matrix of viper scores with samples in columns and regulators
#'  across the rows}
#' \item{cindy}{list of regulators, each with a set of modulators and p values
#' representing their CINDY inferred association}
#' \item{preppi}{list of regulators, each with a set of potential binding
#' partners and PREPPi inferred p values for probability of binding}
#' \item{rawsnp}{matrix of samples and genes with potential mutations.
#' 0 for no mutation, 1 for presence of some non-silent mutation}
#' \item{rawcnv}{matrix of samples and genes with copy number variant scores
#' calculated by GISTIC}
#' \item{fusions}{matrix of samples and fusion genes.
#' 0 for no fusion, 1 for presence of fusion}
#' \item{fCNV}{vector of gene names that have functional CNV mutations}
#' }
"example.gbm.data"