\name{jacquard}
\alias{jacquard}
\title{The Jacquard coefficient as a measure of overlap in terms of genes in 2...}
\usage{jacquard(a, b)
}
\description{The Jacquard coefficient as a measure of overlap in terms of genes in 2
gene sets.}
\details{Motivated by the leading edge analysis of Broad's GSEA.}
\value{a decimal value in [0,1] where 0 represents no overlap, and 1.0
indicating 100\% overlap - ie a=b}
\author{Mark Cowley, 2009-04-06}
\arguments{\item{a}{a vector of entities, be they characters or integers}
\item{b}{a vector of entities, be they characters or integers}
}
\examples{jacquard(letters[1:5], letters[3:5])
# [1] 0.6}

