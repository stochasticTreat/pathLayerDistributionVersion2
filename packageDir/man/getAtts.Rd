\name{getAtts}
\alias{getAtts}
\title{getAtts}
\usage{
getAtts(mat, atts, verbose = F)
}
\arguments{
  \item{mat:}{the matrix to be searched through; contains
  at least 2 columns with column to being the attributes
  column.}

  \item{atts:}{a vector of strings; the attributes to be
  matched. perl compatible regular expressions are
  implemented for string matching}
}
\value{
matrix, the subset of rows matching the strings provided in
atts.
}
\description{
Gets all rows in the matrix, mat, with a the attribute
strings(s) in atts
}

