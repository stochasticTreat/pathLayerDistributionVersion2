\name{checkFileCopyDefault}
\alias{checkFileCopyDefault}
\title{checkFileCopyDefault}
\usage{
checkFileCopyDefault(fname)
}
\arguments{
  \item{fname}{The path of the file to be checked.}
}
\value{
string, the file name.
}
\description{
Checks if a file exists. If it doesn't, function checks for
a default in the package directory's extdata/ folder and
copies the default to the location described by the file
name passed.
}

