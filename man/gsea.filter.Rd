\name{gsea.filter}
\alias{gsea.filter}
\title{Take a GSEA list, and filter out genesets that fail various criteria.}
\usage{gsea.filter(x, N, NES, FDR, P, FWER, direction, genesets)
}
\description{Take a GSEA list, and filter out genesets that fail various criteria.}
\details{The criteria are additive.}
\value{a GSEA list with probably fewer genesets in the top table, the
leading edge list, and the edb object [if present]. also, a
"filter.settings" list is added which details which settings were used}
\author{Mark Cowley, 2009-09-03}
\arguments{\item{x}{a gsea list}
\item{N}{keep the top 'N' genesets, eg 50}
\item{NES}{keep the genesets with |NES| > |NES threshold|.}
\item{FDR}{keep the genesets with NOM.p.val, FDR.q.val or FWER.p.val <
threshold}
\item{P}{keep the genesets with NOM.p.val, FDR.q.val or FWER.p.val <
threshold}
\item{FWER}{keep the genesets with NOM.p.val, FDR.q.val or FWER.p.val <
threshold}
\item{direction}{"up", "down", or NULL or "either"}
\item{genesets}{a vector of geneset ID's such that only these genesets will
be retained.}
}

