\name{getPaths}
\alias{getPaths}
\title{Get set of cellular pathways}
\usage{
getPaths(path_file = NULL,
  referenceFileName = "./reference_data/paths/pathMetaData.txt",
  symtab = NULL, verbose = T)
}
\arguments{
  \item{path_file}{String giving path to file from which
  pathways should be loaded. Alternatively, a \code{Study}
  object can be provided here and pathways will
  automatically be retreived from it.}

  \item{referenceFileName}{Optional. String giving path to
  file where pathway meta data (date loaded, path
  repository name, etc..) can be found.}

  \item{symtab}{A gene identifier lookup table. If not
  provided, one will be automatically obtained (may take
  longer)}

  \item{verbose}{Logical flag indicating if interactive
  user prompts and additional information should be
  displayed.}
}
\value{
\code{Path_Detail} object.
}
\description{
Get \code{Path_Detail} object, housing gene sets for
cellular pathways.
}
\details{
Pathways are stored in the local working directory in the
folder \code{./reference_data/paths/}. Pathways must be
imported into this package so that meta-data and import
records can be maintained. Path meta data is stored in the
file \code{./reference_data/paths/pathMetaData.txt} . Any
time a path repository is imported, a record of its import
is added to the path meta data file. When loading a path
repository for use in a study, this program checks the path
meta data file, for available path repositories and for
path repository file locations.
}
\examples{
\dontrun{
#select and/or load pathways interactively
path_detail = getPaths()
}

#Get pathways from a study object:
study=getTestStudyObject()
path_detail=getPaths(study)

#select pathways from their file name:
path_detail = getPaths(path_file="./reference_data/paths/Reactome.2014.04.06.12.52.27.txt",
											verbose=FALSE)
}

