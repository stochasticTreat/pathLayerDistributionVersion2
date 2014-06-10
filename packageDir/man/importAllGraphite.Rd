\name{importAllGraphite}
\alias{importAllGraphite}
\title{Import all pathway repositories from the graphite package to this package.}
\usage{
importAllGraphite(repositories = list(reactome = reactome, spike = spike, nci
  = nci, kegg = kegg, biocarta = biocarta), symtab = NULL, verbose = F)
}
\arguments{
  \item{repositories}{A \code{list} object, each slot
  containing a pathway repository from graphite that is to
  be imported}

  \item{symtab}{A HUGO symbol look up table, as provided by
  the package HGNChelper or from genenames.org.}

  \item{verbose}{A flag indicating if this function should
  be run in an interactive mode. This allows individual
  pathway repositories to be selected, insted of
  automatically selecting all pathways.}
}
\description{
Gives options to import pathway repositories from Reactome,
Spike, NCI, KEGG and Biocarta.
}

