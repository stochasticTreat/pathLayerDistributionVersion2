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

