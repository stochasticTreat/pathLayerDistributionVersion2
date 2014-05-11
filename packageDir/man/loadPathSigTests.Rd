\name{loadPathSigTests}
\alias{loadPathSigTests}
\title{loadPathSigTests}
\usage{
loadPathSigTests(sigTestFile = "./hypergeometricPathAnalysis.R")
}
\arguments{
  \item{sigTestFile}{An R script file containing one or
  more functions which implement the pathSigTest interface
  and a line defining a named list of variable name
  'pathAnalysisFunctions' with names describing the path
  analysis functions and values set as the respective path
  analysis functions.}
}
\value{
A named list containing the pathway analysis functions to
be applied.
}
\description{
loadPathSigTests
}

