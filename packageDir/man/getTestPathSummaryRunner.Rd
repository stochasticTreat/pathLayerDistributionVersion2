\name{getTestPathSummaryRunner}
\alias{getTestPathSummaryRunner}
\title{Returns an example of the PathSummaryRunner reference class object used by}
\usage{
getTestPathSummaryRunner(pths = NULL)
}
\arguments{
  \item{pths}{A \code{Path_Detail} object. If left null,
  (the default), the default reactome pathways will be
  automatically provided.}
}
\value{
An instance of the \code{PathSummaryRunner} class, loaded
with data indicating enrichment of the "Abacavir
metabolism" Reactome pathway.
}
\description{
Returns an object of the PathSummaryRunner reference class
used by summaryTable.
}

