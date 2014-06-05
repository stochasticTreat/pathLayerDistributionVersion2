\name{getDefaultPaths}
\alias{getDefaultPaths}
\title{Retreive the default set of cellular pathways.}
\usage{

  getDefaultPaths(path_file = system.file("extdata/Reactome.2014.04.06.12.52.27.txt",
  package = "packageDir"))
}
\arguments{
  \item{path_file}{string giving the relative or absolute
  file path to an accessible set of cellular pathways. Set
  of pathways must first have been recorded in the
  reference_data/paths/pathMetaData.txt file.}
}
\value{
A \code{Path_Detail} object.
}
\description{
Retreives the default set of cellular pathways. These are
the gene sets for all Reactome pathways as downloaded on
4/6/14 .
}
\examples{
cellularPathwayGeneSets = getDefaultPaths()
}

