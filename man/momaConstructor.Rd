% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moma.r
\name{MomaConstructor}
\alias{MomaConstructor}
\title{MOMA Constructor Function}
\usage{
MomaConstructor(
  x,
  pathways,
  gene.blacklist = NA_character_,
  output.folder = NA_character_,
  gene.loc.mapping = gene.map,
  viperAssay = "viper",
  mutMat = "mut",
  cnvMat = "cnv",
  fusionMat = "fusion"
)
}
\arguments{
\item{x}{A MultiAssayExerperiment object or list object with the following assays:
(note: by default assays must have these exact names. Otherwise they can be changed
using the viperAssay, mutMat, cnvMat and fusionMat parameters.)
\describe{
\item{viper}{VIPER protein activity matrix with samples as columns 
and rows as protein IDs}
\item{mut}{An indicator matrix (0/1) of mutation events with samples as 
columns and genes as rows}
\item{cnv}{A matrix of CNV scores (typically SNP6 array scores from TCGA) 
with samples as columns and genes as rows}
\item{fusion}{An indicator matrix (0/1) of fusion events with samples as 
columns and genes as rows} }}

\item{pathways}{A named list of lists. Each named list represents 
interactions between proteins (keys) and their associated partners}

\item{gene.blacklist}{A vector of genes to exclude from the analysis}

\item{output.folder}{Location to store output and intermediate results}

\item{gene.loc.mapping}{A data.frame of band locations and Entrez IDs}

\item{viperAssay}{name associated with the viper assay in the assay object}

\item{mutMat}{name associated with the mutation matrix in the assay object}

\item{cnvMat}{name associated with the cnv matrix in the assay object}

\item{fusionMat}{name associated with the fusion matrix in the assay object}
}
\value{
an instance of class Moma
}
\description{
Create MOMA Object from either a MultiAssayExperiment object or a list of
assays.
See vignette for more information on how to set up and run the MOMA object
}
\examples{
momaObj <- MomaConstructor(example.gbm.mae, gbm.pathways)
}
