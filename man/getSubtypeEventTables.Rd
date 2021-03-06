% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.saturation.plots.r
\name{getSubtypeEventTables}
\alias{getSubtypeEventTables}
\title{Helper function to get subtype specific events}
\usage{
getSubtypeEventTables(saturation.data, sample.clustering, checkpoints)
}
\arguments{
\item{saturation.data}{: genomic saturation object from MOMA. List indexed by 
cluster then sample then regulator with the number of events associated with 
each additional regulator}

\item{sample.clustering}{: clustering vector with sample names and 
cluster designations}

\item{checkpoints}{: from momaObj}
}
\value{
a table that has counts of how many times a particular event 
happens in a cluster
}
\description{
Helper function to get subtype specific events
}
\keyword{internal}
