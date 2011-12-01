\docType{methods}
\name{jacquard}
\alias{jacquard}
\alias{jacquard,character,character-method}
\alias{jacquard,GeneSet,GeneSet-method}
\alias{jacquard,numeric,numeric-method}
\title{The Jacquard coefficient as a measure of overlap in terms of genes in 2
gene sets.}
\arguments{
  \item{a}{a vector of entities, be they \code{character}s,
  \code{integer}s or \code{GeneSet}}

  \item{b}{a vector of entities, be they \code{character}s,
  \code{integer}s or \code{GeneSet}}
}
\value{
  a decimal value in [0,1] where 0 represents no overlap,
  and 1.0 indicating 100\% overlap - ie a=b
}
\description{
  Motivated by the leading edge analysis of Broad's GSEA.
}
\examples{
require(GSEABase)
jacquard(letters[1:5], letters[3:5])
jacquard(1:5, 3:5)
}
\author{
  Mark Cowley, 2009-04-06
}

