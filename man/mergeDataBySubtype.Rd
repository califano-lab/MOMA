% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.saturation.plots.r
\name{mergeDataBySubtype}
\alias{mergeDataBySubtype}
\title{Create data frame from coverage data, including number of total events 
'covered' and unique events}
\usage{
mergeDataBySubtype(genomic.saturation, sample.clustering, topN = 100)
}
\arguments{
\item{genomic.saturation}{: data from genomic saturation function}

\item{sample.clustering}{: clustering vector with sample names and 
cluster designations}

\item{topN}{: number of regulators to look through. default is 100}
}
\value{
dataframe with coverage data for genomic events
}
\description{
Create data frame from coverage data, including number of total events 
'covered' and unique events
}
\keyword{internal}
