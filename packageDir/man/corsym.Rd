\name{corsym}
\alias{corsym}
\title{Correct unofficial gene symbols.}
\usage{
corsym(symbol_set, symref = NULL, verbose = T, col2 = "Chrom",
  correctionsfile = "./reference_data/gene_symbol_corrections_list.txt")
}
\arguments{
  \item{symbol_set}{the set of gene identifiers to be
  corrected}

  \item{symref}{The gene identifier lookup table containing
  a set of official gene symbols. If this contains two
  columns, it should be in the format of the hgnc.table
  data set provided by the CRAN package, HGNChelper. Else,
  this should be a whole data table of gene symbol look ups
  as provided by genenames.org}

  \item{verbose}{If set to TRUE, this flag will cause the
  gene symbol corrections to be conducted automatically.}

  \item{col2}{If two columns of data are provided as the
  symbol_set input to corsym, col2 should provide the name
  of a second column containing additional data bout the
  gene identifier to be checked, such as a chromosomal
  location.}

  \item{correctionsfile}{character string giving the file
  path to the a gene symbol corrections file. This file
  should contain two columns: old_symbol and new_symbol
  containing the errant and the correct symbols,
  respectively. This is the primary set of symbols to be
  used in coordinating gene symbol corrections between the
  multiple arms of the study.}
}
\value{
The character vector of gene identifiers provided as an
input, with any possible corrections made.
}
\description{
Function to correct and unify gene identifiers across a
study arms.
}
\details{
This function allows a user to interactively correct gene
identifiers so they can be coordinated between study arms.
All corrections can optionally be recorded to the
corrections file so that at latter time they can be run
automatically, without user interaction, and so that they
can be automatically re-used in correcting gene symbols
from other study arms.
}
\examples{
toCheck = c("p53", "FLT3", "ASM1", "ASML3B","ARF","APOBEC3C")
corrected = corsym(symbol_set=toCheck) #using HGNChelper and previously made corrections

symref  = getHugoSymbols()

\dontrun{
symref  = getHugoSymbols(curhugofname="./reference_data/current_hugo_table_slim.txt")
corrected = corsym(symbol_set=toCheck, symref=symref, verbose=T)
}
}

