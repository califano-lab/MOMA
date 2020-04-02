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