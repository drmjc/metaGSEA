\name{gsea.leadingedge.distance}
\alias{gsea.leadingedge.distance}
\title{Determine the distance between various leading edge genes, in terms of the...}
\usage{gsea.leadingedge.distance(x)
}
\description{Determine the distance between various leading edge genes, in terms of the
jacquard.}
\value{an object of class \sQuote{dist} where each value is \code{1-Jacquard(a, b)}
where \code{a} and \code{b} are the leading edges from two genesets.}
\author{Mark Cowley, 2009-04-06}
\seealso{\code{\link{dist}}, \code{\link{dist.gsea}}}
\arguments{\item{x}{a list of vectors of gene symbols}
}

