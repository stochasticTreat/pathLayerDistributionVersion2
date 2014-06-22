\name{getReactomeBiopax}
\alias{getReactomeBiopax}
\title{Download biopax files for reactome pathways.}
\usage{
getReactomeBiopax(study, pathNames, verbose = T)
}
\arguments{
  \item{study}{A \code{Study} object.}

  \item{pathNames}{The names of the pathways the program
  should attempt to download.}

  \item{verbose}{\code{logical} flag indicating if user
  should be prompted before attempted download of biopax
  pathways.}
}
\value{
Character vector containing the names of any pathways that
could not be downloaded.
}
\description{
Attempts to download biopax files for reactome pathways.
First, pathway's individual database identifiers are
obtained from bioMART, then restful calls are made to
Reactome's restful interface. A record of successfully
downloaded pathways is made in the file
"./reference_data/paths/biopax/record_of_biopax_pathways.txt"
so that pathways will not be repeatedly downloaded.
}

