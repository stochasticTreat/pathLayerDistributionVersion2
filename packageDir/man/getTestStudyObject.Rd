\name{getTestStudyObject}
\alias{getTestStudyObject}
\title{Get a Study object containing moc data.}
\usage{
getTestStudyObject(noSettings = FALSE, noResults = FALSE, noPaths = FALSE,
  noArms = FALSE)
}
\arguments{
  \item{noSettings}{\code{logical} flag indicating if
  analysis setting should not be automatically added to
  \code{Study} object.}

  \item{noResults}{\code{logical} flag indicating if
  analysis results should not be automatically added to
  \code{Study} object.}

  \item{noPaths}{\code{logical} flag indicating if cellular
  pathway repository should not be automatically added to
  \code{Study} object.}

  \item{noArms}{\code{logical} flag indicating if data
  processing arms should not be loaded with \code{Study}
  object.}
}
\value{
\code{Study} object
}
\description{
Returns a \code{Study} object loaded with moc somatic
mutation data which indicates the Reactome pathway,
Abacavir metabolism, is enriched in somatic mutation.
}

