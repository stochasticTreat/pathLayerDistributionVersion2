\name{getPublicationRefs}
\alias{getPublicationRefs}
\title{getPublicationRefs()}
\usage{
getPublicationRefs(df, tdf)
}
\arguments{
  \item{df}{the main biopax data frame}

  \item{tdf}{a subsegment of the biopax data frame from
  which citations and citation dates are to be found. Must
  have columns "property", "property_value" and "id".}
}
\value{
data frame with columns "id" "citation" "date"
}
\description{
getPublicationRefs()
}

