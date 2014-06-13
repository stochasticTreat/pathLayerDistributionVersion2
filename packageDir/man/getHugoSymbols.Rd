\name{getHugoSymbols}
\alias{getHugoSymbols}
\title{getHugoSymbols}
\usage{
getHugoSymbols(paths_detail = NULL, curhugofname = NULL, verbose = F)
}
\arguments{
  \item{paths_detail}{A \code{Path_Detail} or a
  \code{Study} object. If a \code{Path_Detail} object is
  provided, the symbols will be extracted from its
  \code{symtable} slot. If a \code{Study} object is
  provided for this argument, the symbol lookup table will
  be exracted from the \code{Path_Detail} it contains.}

  \item{curhugofname}{If interactive gene symbol correction
  is to be used, this argument should be the file path to
  the HUGO table as downloaded from genenames.org.}

  \item{verbose}{Controlls if symbol corrections are to be
  interactive (if yes, curhugofname file must be supplied
  as it contains critical information, such as gene symbol
  status, past identifiers and synonyms)}
}
\value{
Table of symbols: either a two column data.frame, the
hgnc.table provided by HGNChelper, or a data frame as
provided by genenames.org. Either of which have a column
titled Approved.Symbol which contains official, approved
symbols.
}
\description{
obtains HGNC/HUGO symbol look up table containing official
set of HUGO symbols.
}
\examples{
#get the default HUGO lookup table, hgnc.table from HGNChelper.
hsyms = getHugoSymbols()
#attempt to download and or use a full hugo lookup table from genenames.org
hsyms = getHugoSymbols(curhugofname="./reference_data/current_hugo_table_slim.txt")
}

