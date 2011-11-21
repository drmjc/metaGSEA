\name{GMThyperGTest}
\alias{GMThyperGTest}
\title{Perform a Hypergeometric test, upon a GMT file.
Much like \code{GOhyperGTest} from \code{topGO}, this function computes
the overlaps between a vector
of gene id's of interest, and a collection of gene sets, stored within a gmt file}
\usage{
  GMThyperGTest(genes, gmt, adjust.method = "none")
}
\arguments{
  \item{genes}{a vector of gene symbols}

  \item{gmt}{a gmt object. see
  \code{\link{import.gsea.gmt}}}

  \item{adjust.method}{the multiple hypothesis correction
  value. see \code{\link{p.adjust}}.
  default=\dQuote{none}}
}
\value{
  a \code{data.frame} containing these columns:
  \dQuote{GeneSet}, \dQuote{x}, \dQuote{m}, \dQuote{n},
  \dQuote{k}, \dQuote{Pval}, \dQuote{FDR}, where:\cr
  \dQuote{x} = the number of white balls drawn without
  replacement from an urn which contains both black and
  white balls.\cr \dQuote{m} = the number of white balls in
  the urn.\cr \dQuote{n} = the number of black balls in the
  urn.\cr \dQuote{k} = the number of balls drawn from the
  urn.\cr
}
\description{
  The type of identifier used doesn't matter much, as long
  as the same identifiers are used in genes and gmt, eg
  Gene Symbols.
}
\author{
  Mark Cowley, 2009-01-15
}
\seealso{
  \code{\link[stats]{phyper}}
}

