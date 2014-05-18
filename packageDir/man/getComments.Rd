\name{getComments}
\alias{getComments}
\title{getComments()}
\usage{
getComments(tdf)
}
\arguments{
  \item{tdf}{data.frame, the section of the biopax data
  frame to pull comments about nodes from. Must have
  columns "property", "property_value" and "id"}
}
\value{
data frame with two columns "id" and "comment" where
multiple comments are sepparated by pipes (ie, "|")
}
\description{
Gets and pastes together comments for all elements in the
data.frame provided by rBiopaxParser's readBiopax function.
}

