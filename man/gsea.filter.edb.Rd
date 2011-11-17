\name{gsea.filter.edb}
\alias{gsea.filter.edb}
\title{Filter a GSEA edb file, restricting the entries to restricted set of
genesets.}
\usage{
  gsea.filter.edb(x, genesets)
}
\arguments{
  \item{x}{a gsea list (not just the edb file)}

  \item{genesets}{a vector of geneset names. Minimal error
  checking done, so make sure they exist in your edb file!}
}
\value{
  a gsea list with the edb object restricted to only those
  that were specified by \code{genesets}
}
\description{
  Filter a GSEA edb file, restricting the entries to
  restricted set of genesets.
}
\author{
  Mark Cowley, 2009-10-13
}

