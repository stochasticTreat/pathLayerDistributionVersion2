\name{prepPathListForSave}
\alias{prepPathListForSave}
\title{prepPathListForSave()}
\usage{
prepPathListForSave(gset, path_source = "Reactome pathway")
}
\arguments{
  \item{gset}{The list object of paths}

  \item{path_source}{A string desribing the source of the
  pathways, to be pasted as the second element in each line
  of the output in the manner as "Reactome pathway" is seen
  here: Abacavir metabolism Reactome pathway UniProt:P49902
  UniProt:Q16774}
}
\value{
one column matix, with each row containing the
tab-delimited set: path name, path source and path members
}
\description{
preps paths in gene-set format (list with names=path names
and values = vector of gene names) for saving as a tab
delimited file
}
\examples{
library(graphite)
plistsAsDataFrame = prepPathListForSave(gset=reactome, path_source="Reactome pathway")
}

