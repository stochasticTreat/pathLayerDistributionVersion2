\name{addCorrections}
\alias{addCorrections}
\title{Function to add corrections to the corrections file. Internal, used by corsym to coordinate symbol corrections.}
\usage{
addCorrections(new_corrections,
  correctionsfile = "./reference_data/gene_symbol_corrections_list.txt",
  correction_set = NULL)
}
\arguments{
  \item{new_corrections}{a two column matrix; column 1 =
  old, incorrect symbols, column 2 = new, corrected
  symbols}

  \item{correctionsfile}{the file name of the corrections
  file}

  \item{correction_set}{the original contents of the
  corrections file before new corrections were made}
}
\description{
Adds corrections to the corrections file
}

