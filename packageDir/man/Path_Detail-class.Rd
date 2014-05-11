\docType{class}
\name{Path_Detail-class}
\alias{Path_Detail-class}
\title{Path_Detail class}
\description{
Contains all information relevent to cellular pathways used
in pathway analysis.
}
\section{Slots}{

\describe{
\item{\code{name}}{character. The name of the pathway repository used}

\item{\code{file}}{character. The file from which the pathway repostory was loaded}

\item{\code{date}}{character. The date the pathway repository was loaded}

\item{\code{info}}{character. Any additional information reguarding the pathway repository.}

\item{\code{source}}{character. The source of the pathway repository. (ex: Reactome, via the graphite bioconductor package )}

\item{\code{graphite.}}{Depricated. Only here for legacy purposes.}

\item{\code{gene_overlap_counts}}{Number of pathways each gene belongs to.}

\item{\code{full_path_length}}{The number of genes annotated to each pathway.}

\item{\code{symtable}}{The official, approved set of symbols or gene identifiers and any cross references.}

\item{\code{paths}}{matrix. The actual storage of pathways in bipartate graph format, with row names as the path names, column names as gene identifiers and values logical representing gene membership in pathawys.}

\item{\code{symbol_type.}}{character. The type of gene identiers used, for example HUGO or Uniprot.}
}}

