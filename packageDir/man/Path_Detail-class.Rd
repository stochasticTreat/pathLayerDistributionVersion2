\docType{class}
\name{Path_Detail-class}
\alias{Path_Detail-class}
\title{Path_Detail class}
\description{
Contains all information relevent to cellular pathways used
in pathway analysis.
}
\section{Fields}{
  \describe{ \item{\code{name}:}{Object of class
  \code{"character"} The name of the pathway repository
  used. } \item{\code{file}:}{Object of class
  \code{"character"} The file from which the pathway
  repostory was loaded.} \item{\code{date}:}{Object of
  class \code{"character"} The date the pathway repository
  was loaded.} \item{\code{info}:}{Object of class
  \code{"character"} Any additional information reguarding
  the pathway repository.} \item{\code{source}:}{Object of
  class \code{"character"} The source of the pathway
  repository. For example, Reactome, via the graphite
  bioconductor package )}
  \item{\code{graphite}:}{Depricated. Only here for legacy
  purposes.} \item{\code{gene_overlap_counts}:}{Object of
  class \code{"matrix"} One column, with row names given as
  gene names. Values are the number of pathways each gene
  belongs to.} \item{\code{full_path_length}:}{Object of
  class \code{"matrix"} One column, pathway names given as
  row names. Values given as the number of genes annotated
  to each pathway.} \item{\code{symtable}:}{Object of class
  \code{"data.frame"} The official, approved set of symbols
  or gene identifiers and any cross references.}
  \item{\code{paths}:}{Object of class \code{"matrix"} The
  actual storage of pathways in bipartate graph format,
  with row names as the path names, column names as gene
  identifiers and values logical representing gene
  membership in pathawys.}
  \item{\code{symbol_type}:}{Object of class
  \code{"character"} The type of gene identiers used, for
  example HUGO or Uniprot.} }
}

