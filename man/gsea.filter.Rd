\name{gsea.filter}
\alias{gsea.filter}
\title{Take a GSEA list, and filter out genesets that fail various criteria.}
\usage{
  gsea.filter(x, N = NULL, NES = NULL, FDR = NULL,
    P = NULL, FWER = NULL, direction = NULL,
    genesets = NULL)
}
\arguments{
  \item{x}{a gsea list}

  \item{N}{keep the top 'N' genesets, eg 50}

  \item{NES}{keep the genesets with |NES| > |NES
  threshold|.}

  \item{FDR}{keep the genesets with NOM.p.val, FDR.q.val or
  FWER.p.val < threshold}

  \item{P}{keep the genesets with NOM.p.val, FDR.q.val or
  FWER.p.val < threshold}

  \item{FWER}{keep the genesets with NOM.p.val, FDR.q.val
  or FWER.p.val < threshold}

  \item{direction}{"up", "down", or NULL or "either"}

  \item{genesets}{a vector of geneset ID's such that only
  these genesets will be retained.}
}
\value{
  a GSEA list with probably fewer genesets in the top
  table, the leading edge list, and the edb object [if
  present]. also, a \code{filter.settings} list is added
  which details which settings were used
}
\description{
  The criteria are additive.
}
\author{
  Mark Cowley, 2009-09-03
}

