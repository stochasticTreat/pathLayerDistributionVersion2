\name{downloadHugoLookupTable}
\alias{downloadHugoLookupTable}
\title{Obtain full HUGO gene symbol look up table.}
\usage{

  downloadHugoLookupTable(curhugofname = "./reference_data/current_hugo_table_slim.txt",
  autoDownLoad = F)
}
\arguments{
  \item{curhugofname}{\code{character string} The location
  the symbol lookup table should be saved. Default is the
  local directory at
  ./reference_data/current_hugo_table_slim.txt. For this
  program to use the full HUGO table, the file needs to be
  left here.}

  \item{autoDownLoad}{\code{logical} A flag indicating if
  the program should automatically attempt to download the
  table (T) or if it should just display information on how
  to manually download the table (F).}
}
\description{
The function attempts to download the HUGO lookup table
(includes previous symbols, synonyms and chromosome
corrdinates) from the HGNC website genenames.org.
}

