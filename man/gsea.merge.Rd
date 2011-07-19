\name{gsea.merge}
\alias{gsea.merge}
\title{Merge two GSEA objects.}
\usage{gsea.merge(gsea1, gsea2, prefix, suffix, warn=TRUE)
}
\description{Merge two GSEA objects.}
\details{You would want to do this if you compare two rnk files to 1 gmt file (eg
c2_all), and want to compare the significant genesets between these runs.
It's no good knowing that geneset X is in the top 5 in both runs if they
contain markedly different genes in their leading edges.}
\author{Mark Cowley, 2009-10-14}
\arguments{\item{gsea1}{GSEA object to be merged.}
\item{gsea2}{GSEA object to be merged.}
\item{prefix}{character(2) containing optional prefix; or \code{NULL}}
\item{suffix}{character(2) containing optional suffix; or \code{NULL}}
\item{warn}{logical: produce a warning if \code{prefix=NULL && suffix=NULL}}
}

