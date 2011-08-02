\name{fisherMethod.GSEAlist}
\alias{fisherMethod.GSEAlist}
\title{Apply Fisher's method to a GSEA list...}
\usage{fisherMethod.GSEAlist(gsea.list, values=c("NOM.p.val", "FDR.q.val"), direction=c("either",
    "same"))
}
\description{Apply Fisher's method to a GSEA list}
\details{Fisher's method is a form of meta-analysis, where a single p-value is
computed for a set of input p-values. This uses the p-values within the
top table objects.
The GSEA results must have been run against the same GMT file.}
\value{a data.frame with ID, {S, p.value, q.value, rank} all produced by
fisherMethod, and the p-values from each of the N GSEA results.}
\author{Mark Cowley, 2010-07-22}
\arguments{\item{gsea.list}{a list of GSEA objects}
\item{values}{Which statistic to use, the unadjusted P-values (\dQuote{NOM.p.val}), 
or the FDR q values (\dQuote{FDR.q.val})?}
\item{direction}{either: ignore the direction of the gene's change.}
}

