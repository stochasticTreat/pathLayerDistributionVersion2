\name{makeDrugSelectionWorksheet}
\alias{makeDrugSelectionWorksheet}
\title{Uses genomic data from the provided Study object to produce a table of pertinent drug-gene associations.}
\usage{
makeDrugSelectionWorksheet(STUDY, pathsToSearch = NULL)
}
\arguments{
  \item{STUDY}{A \code{Study} object}

  \item{pathsToSearch}{A \code{vector} of pathway names.
  The set of cellular pathways that is to be searched for
  drug targets.}
}
\value{
Spread sheet relating genes in paths to drugs and clinical
trials.
}
\description{
Uses genomic data from the provided Study object to produce
a table of pertinent drug-gene associations.
}

