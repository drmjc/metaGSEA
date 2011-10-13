\name{jacquard-methods}
\alias{jacquard-methods}
\alias{jacquard}
\title{The Jacquard coefficient as a measure of overlap in terms of genes in 2...}
\description{The Jacquard coefficient as a measure of overlap in terms of genes in 2
gene sets.}
\details{\code{jacquard-methods}: Motivated by the leading edge analysis of Broad's GSEA.

}
\value{\code{jacquard-methods}: a decimal value in [0,1] where 0 represents no overlap, and 1.0
indicating 100\% overlap - ie a=b

}
\author{Mark Cowley, 2009-04-06}
\docType{methods}
\arguments{\item{a}{a vector of entities, be they \code{character}s, \code{integer}s or \code{GeneSet}}
\item{b}{a vector of entities, be they \code{character}s, \code{integer}s or \code{GeneSet}}
}
\examples{jacquard(letters[1:5], letters[3:5])
jacquard(1:5, 3:5)}
\alias{jacquard,character,character-method}
\alias{jacquard}
\usage{\S4method{jacquard}{character,character}(a, b)
\S4method{jacquard}{numeric,numeric}(a, b)
\S4method{jacquard}{GeneSet,GeneSet}(a, b)
}
\alias{jacquard,character,character-method}
\alias{jacquard,numeric,numeric-method}
\alias{jacquard}
\alias{jacquard,numeric,numeric-method}
\alias{jacquard,GeneSet,GeneSet-method}
\alias{jacquard}
\alias{jacquard,GeneSet,GeneSet-method}

