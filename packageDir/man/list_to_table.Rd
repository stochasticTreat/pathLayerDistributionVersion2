\name{list_to_table}
\alias{list_to_table}
\title{Convert list of vectors to bipartate graph format.}
\usage{
list_to_table(pth, uni = NULL)
}
\arguments{
  \item{pth}{A list object, with each slot containing a
  vector.}

  \item{uni}{The full, unique set of elements in all the
  vectors in the list passed in the argument, pth. If not
  provided this will be automatically recalculated.}
}
\value{
\code{matrix} object with logical values representing a
bipartate graph; rows are set as the slot names of the
input list. Columns are set as all the unique elements in
the vectors the list contains.
}
\description{
Used to restructure list of vectors to bipartate graph
format.
}

