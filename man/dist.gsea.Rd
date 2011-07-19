\name{dist.gsea}
\alias{dist.gsea}
\title{Determine the distance between various leading edge genes, in terms of the...}
\usage{dist.gsea(x)
}
\description{Determine the distance between various leading edge genes, in terms of the
jacquard.}
\details{This is a convenience wrapper around \code{\link{gsea.leadingedge.distance}}}
\value{an object of class 'dist' where each value is 1-Jacquard(a, b)
where a and b are the leading edges from two genesets.}
\author{Mark Cowley, 2009-04-06}
\seealso{\code{\link{dist}}, \code{\link{gsea.leadingedge.distance}}}
\arguments{\item{x}{a list of vectors of gene symbols}
}

