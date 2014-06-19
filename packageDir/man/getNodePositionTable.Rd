\name{getNodePositionTable}
\alias{getNodePositionTable}
\title{Get position of nodes in cytoscape diagram. .}
\usage{
getNodePositionTable(nodeTable, attval, attname, w)
}
\arguments{
  \item{nodeTable}{The table of nodes}

  \item{attval}{Used to get positions of a subset of nodes.
  Ex: if the attname is cellularLocation, and the attval is
  cell_membrane, the nodes with the cellularLocation
  attribute set to cell_membrane would be returned.}

  \item{attname}{The name of the attribute to filter the
  nodes by}

  \item{w}{The cytoscape window connection object.}
}
\value{
data table with three columns: node id, x position, y
postion
}
\description{
Gets a table with the positions of all the nodes.
}

