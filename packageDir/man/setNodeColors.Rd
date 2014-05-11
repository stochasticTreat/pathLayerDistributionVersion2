\name{setNodeColors}
\alias{setNodeColors}
\title{setNodeColors()}
\usage{
setNodeColors(w, study, node_color_choices, targeted, active = NULL,
  default.color = NULL)
}
\arguments{
  \item{w:}{the cytoscape connection}

  \item{node_color_choices:}{a vector, length 3, containing
  the hex color codes for the targeted, the unactive and
  the active nodes}

  \item{targeted:}{the vector of nodes considered to be
  targeted}

  \item{active:}{the vector of nodes considered to be
  active}

  \item{default.color:}{the hex color used as the default
  color}
}
\description{
setNodeColors()
}

