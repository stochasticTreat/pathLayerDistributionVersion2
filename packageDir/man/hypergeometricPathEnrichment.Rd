\name{hypergeometricPathEnrichment}
\alias{hypergeometricPathEnrichment}
\title{hypergeometric Path Enrichment analysis function for summaryTable interface.}
\usage{
hypergeometricPathEnrichment(pathSigRunner, paths_detail)
}
\arguments{
  \item{paths_detail}{An object of the Path_Detail class}

  \item{pathSigRunner}{An object of the PathSummaryRunner
  reference class}
}
\value{
data.frame with two columns, one with hypergeometric
p-values for each path, the other with B.H./FDR adjusted
p-value.
}
\description{
Provides hypergeometric test for analysis of pathway
enrichment. Implements interface from the summaryTable
function to path significance tests.
}
\note{
Implements the path test interface with one input, the path
summary list, and the output, a data.frame with pathway
names as row name(s) and named column(s) giving the result
of the significance test.
}

