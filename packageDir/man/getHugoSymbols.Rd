\name{getHugoSymbols}
\alias{getHugoSymbols}
\title{getHugoSymbols}
\usage{
getHugoSymbols(paths_detail = NULL, curhugofname = NULL, verbose = T)
}
\arguments{
  \item{paths_detail}{If a \code{Path_Detail} object is
  provided, the symbols will be extracted from its
  \code{symtable} slot}

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

