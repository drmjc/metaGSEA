\name{chip2description}
\alias{chip2description}
\title{Convert a GSEA chip file or object into a vector of descriptions.}
\usage{
  chip2description(chip.file = NULL, chip = NULL,
    sep = ": ", probes = NULL, genoloc = NULL)
}
\arguments{
  \item{chip.file}{the path to a chip file; or:}

  \item{chip}{an imported chip object (see
  \code{\link{import.gsea.chip}})}

  \item{sep}{the field separator between the symbol and
  description}

  \item{probes}{an optional vector of probe names to limit
  the result to. Default is \code{NULL} whereby all probes
  will be converted}

  \item{genoloc}{optional vector of genomic locations,
  \code{NA} and \code{NULL} allowed, MUST be in identical
  order to the rows of the chip}
}
\value{
  a character vector of description. If genoloc is
  provided, then eg \dQuote{|@chr6:1950428-1950681|} will
  be appended.
}
\description{
  Convert a GSEA chip file or object into a vector of
  descriptions.
}
\author{
  Mark Cowley
}

