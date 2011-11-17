\name{dist.gsea}
\alias{dist.gsea}
\title{Determine the distance between various leading edge genes, in terms of the
jacquard.}
\usage{
  dist.gsea(x)
}
\arguments{
  \item{x}{a list of vectors of gene symbols}
}
\value{
  an object of class 'dist' where each value is
  1-Jacquard(a, b) where a and b are the leading edges from
  two genesets.
}
\description{
  This is a convenience wrapper around
  \code{\link{gsea.leadingedge.distance}}
}
\author{
  Mark Cowley, 2009-04-06
}
\seealso{
  \code{\link{dist}},
  \code{\link{gsea.leadingedge.distance}}
}

