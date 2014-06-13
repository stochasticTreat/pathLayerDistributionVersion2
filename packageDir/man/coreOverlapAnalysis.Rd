\name{coreOverlapAnalysis}
\alias{coreOverlapAnalysis}
\title{The central function used in an overlap analysis}
\usage{
coreOverlapAnalysis(functionalEnrichmentAnalysis, combinedAberrations, ola,
  coverage = NULL)
}
\arguments{
  \item{functionalEnrichmentAnalysis}{The functional
  enrichment analysis to be used.}

  \item{combinedAberrations}{The aberration analysis to be
  used}

  \item{ola}{The overlap analysis runner object to be
  used.}

  \item{coverage}{If coverage is limited, this will contain
  the set of genes the coverage is limited to.}
}
\value{
An overlap analysis
}
\description{
Examines the actual overlap between functionally
significant and aberrationally signficant pathways.
}

