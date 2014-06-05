\name{setCutoff}
\alias{setCutoff}
\title{Set the cutoff for drug screen score}
\usage{
setCutoff(pgm, s, verbose = T, cutoffScore = NULL)
}
\arguments{
  \item{pgm}{A patient gene matrix: bipartate graph with
  patient ids as columns, genes as rows and cell values the
  drug screen score.}

  \item{s}{A settings list object.}

  \item{verbose}{A logical flag indicating if additional
  information should be displayed.}

  \item{cutoffScore}{Optional, the score, above which genes
  will be considered sensitive. If not provided,
  interactive prompt will be shown to the user.}
}
\value{
\code{list} with slots pgm, a logical patient gene matrix
and \code{s} the settings list object.
}
\description{
Anything above cutoff will be considered to be
drug-sensitive.
}

