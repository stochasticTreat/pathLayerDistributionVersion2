\name{correctByHgncHelper}
\alias{correctByHgncHelper}
\title{Correct gene symbols using the HGNChelper package}
\usage{
correctByHgncHelper(symbol_set, symtab, correctionsfile, verbose = F)
}
\arguments{
  \item{symbol_set}{Vector, matrix or data.frame (if more
  than 2 dimensions, the symbols must be in the first
  column)}

  \item{symtab}{The table of official, orthodox gene
  identifiers.}

  \item{correctionsfile}{A file name for a file of gene
  symbol corrections. This file should be tab delimited and
  contain two named columns: old_symbol and new_symbol. The
  old_symbol column contains incorrect identifiers to be
  corrected, and the new_symbol column contains the
  corrections. There should be no row names.}

  \item{verbose}{A logical flag indicating if interactive
  mode should be run.}
}
\value{
A character vector containing the identifiers for the set
of genes provided in the symbol_set function argument, but
with any applicable corrections made.
}
\description{
Correct gene symbols using the HGNChelper package. Provides
interactive or automatic correction of gene identifiers,
and corrdination of any correction between multiple calls
to this functions.
}

