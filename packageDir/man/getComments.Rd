\name{getComments}
\alias{getComments}
\title{getPublicationRefs()
#}
\usage{
getComments(tdf)
}
\arguments{
  \item{df:}{the main biopax data frame #}

  \item{tdf:}{a subsegment of the biopax data frame from
  which citations and citation dates are to be found. Must
  have columns "property", "property_value" and "id". #}

  \item{tdf:}{data frame, the section of the biopax data
  frame to pull comments about nodes from. Must have
  columns "property", "property_value" and "id"}
}
\value{
data frame with columns "id" "citation" "date"

data frame with two columns "id" and "comment" where
multiple comments are sepparated by pipes (ie, "|")
}
\description{
#
}

