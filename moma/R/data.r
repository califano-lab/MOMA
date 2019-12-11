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
#' \item{gene.loc.mapping}{dataframe with 4 columns, Gene Symbol, Locus ID, 
#' Cytoband location, Entrez ID. Used for mapping amplifications and deletions 
#' between gene and location}
#' \item{clinical}{dataframe of clinical annotations, taken directly from TCGA}
#' \item{fCNV}{vector of gene names that have functional CNV mutations}
#' \item{mutSig}{blacklist of genes to not use in the DIGGIT analysis}
#' }
"gbm.example"

#' MutSig Blacklisted genes
#' 
#' List of genes to not include in the DIGGIT mutation inference because they 
#' have been found to be mutated more often than expected by chance given 
#' background mutation processes.
#' 
#' @format A character vector of Entrez Gene IDs 
#' @source \url{https://software.broadinstitute.org/cancer/cga/mutsig}
"mutSig"

#' Gene Location Mapping
#' 
#' Table used for converting between different forms of gene information.
#' Downloaded from HGNC's custom download portal using the "Approved Symbol", 
#' "NCBI Gene ID", "Chromosome" and "Ensembl Gene ID" curated data options 
#' and only those with "Approved" status.
#' Updated December 2019.
#' 
#' @format A Data frame with 4 columns
#' \describe{
#' \item{Gene.Symbol}{Approved Symbol gene name}
#' \item{Entrez.IDs}{NCBI Gene ID}
#' \item{Cytoband}{Chromosome location}
#' \item{Ensembl}{Ensembl gene ID}
#' }
#'  @source \url{https://www.genenames.org/download/custom/}
"gene.map" 