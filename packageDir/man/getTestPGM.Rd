\name{getTestPGM}
\alias{getTestPGM}
\title{Get an example patient gene matrix}
\usage{
getTestPGM()
}
\value{
\code{matrix} object with column names set as patient IDs,
row names set as HUGO gene symbols, and \code{logical} cell
values indicating if the corresponding gene is active in
the corresponding patient.
}
\description{
Retreives from package data a patient gene matrix with
corresponds with the Reactome pathway, 'Abacavir
metabolism' being completely enriched in active (mutated,
drug sensitive, etc.) genes.
}

