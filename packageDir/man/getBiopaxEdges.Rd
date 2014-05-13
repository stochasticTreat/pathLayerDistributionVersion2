\name{getBiopaxEdges}
\alias{getBiopaxEdges}
\title{getBiopaxEdges()}
\usage{
getBiopaxEdges(df, connectionTypes, classType)
}
\arguments{
  \item{df}{the biopax data frame}

  \item{connectionTypes}{string, the type of
  connection/edge to be searched for}

  \item{classType}{string, the class of nodes be searched
  through}
}
\value{
matrix with columns type, participant and direction: the
type of edge, the participant and the direction or type of
interaction
}
\description{
retreives edge data from the data frame part of a biopax
object (the type obtained from RbiopaxParser)
}

