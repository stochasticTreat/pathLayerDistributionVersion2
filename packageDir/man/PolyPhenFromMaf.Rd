\name{PolyPhenFromMaf}
\alias{PolyPhenFromMaf}
\title{orchestrates creation of polyphen input file from .maf file}
\usage{
PolyPhenFromMaf(mafFname = NULL, outFname = NULL)
}
\arguments{
  \item{mafFname}{character string. The file name and full
  or relative path to the maf input file}

  \item{outFname}{character string. The name of the file
  for the polyphen input file to be saved to.}
}
\value{
A list of diagnostic information, with each slot containing
a specific data point, such as "Rows of mutation data found
in .maf file input"
}
\description{
Allows user to create valid polyphen input file from TCGAs
somatic mutation data storage format, .maf.
}

