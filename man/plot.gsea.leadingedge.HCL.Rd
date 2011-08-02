\name{plot.gsea.leadingedge.HCL}
\alias{plot.gsea.leadingedge.HCL}
\title{GSEA Leading Edge HCL plot.}
\usage{plot.gsea.leadingedge.HCL(x, N, NES, FDR, P, FWER, direction=c("either", "up", "down"),
    main="Leading edge similarities", xlab="", yaxis=c("similarity",
    "distance"), max.distance=1, rename.genesets=FALSE,
    label.clusters=FALSE, cluster.threshold=0.99,
    hclust.method="complete", ...)
}
\description{GSEA Leading Edge HCL plot.}
\details{This function compares genesets in terms of the overlap of the leading edge
genes
within each geneset. This information is presented as a hierarchical
clustering plot,
where distance is measured proportional to the number of genes in common
(see \code{\link{jacquard}}).
The user can either pass in a leading edge list (ie a list where each
element is
a vector of gene symbols), or (as will be more common),
an entire leading edge list (eg 1434 elements for c2_all), in conjunction
with a GSEA
top table (tt), coupled with a threshold, such as the top 50, or those with
FDR < 0.25
via the N, FDR, P, or FWER arguments.}
\note{max.distance\cr
The intuitive outcome of setting \code{max.distance} (eg 0.1)
is that no genesets with a distance of >0.1 (ie a similarity of <90\%)
should remain in the HCL. This is not quite the case. We filter out genesets
that have distance to all other genesets > 0.9 in the 'distance' object
(see \code{\link{gsea.leadingedge.distance}}) which is slightly different to filtering
out nodes in the resulting HCL (which is made upon the distance object).
Net effect:\cr
- You can get singletons in the filtered HCL because all of
the things that the geneset was connected to have been filtered out.\cr
- You can get genesets in a tree of significant genesets, which is itself
below the distance threshold due to a similar explanation to above.}
\value{invisibly returns the hclust object.}
\author{Mark Cowley, 2009-04-06}
\seealso{\code{\link{plot.gsea.leadingedge.HCL}} \code{\link{plclust.gsea}} \code{\link{gsea.leadingedge.distance}}}
\arguments{\item{x}{a GSEA object. see \code{\link{import.gsea}}}
\item{N}{if not \code{NULL}, then choose the top N most significant genesets}
\item{NES}{if not \code{NULL}, then choose those genesets with |NES| > the
threshold}
\item{FDR}{if not \code{NULL}, then choose those genesets with FDR < the
threshold}
\item{P}{if not \code{NULL}, then choose those genesets with nominal P <
the threshold}
\item{FWER}{if not \code{NULL}, then choose those genesets with FWER < the
threshold}
\item{direction}{include just the genesets that are \dQuote{up}, \dQuote{down}, or
\dQuote{either} means top 50 genesets that are either up or down}
\item{main}{additional arguments passed to plot.}
\item{xlab}{additional arguments passed to plot.}
\item{yaxis}{either \dQuote{similarity} or \dQuote{distance}. 
Two closely related genesets have eg 90\% similarity, and a distance of 0.1.}
\item{max.distance}{(experimental) attempt to filter out genesets that are
dissimilar to all other genesets. This is a distance threshold, not a
similarity threshold, so values of 0.9 are a good place to start. See
note.}
\item{rename.genesets}{should the genesets be renamed to include the
rank/fdr/direction}
\item{label.clusters}{if TRUE, then draw red boxes around the clusters. see
cluster.threshold.}
\item{cluster.threshold}{the distance threshold if label.clusters=TRUE.
0.99 is a good default in practice.}
\item{\dots}{additional arguments passed to plot.}
\item{hclust.method}{the hierarchical clustering method. see \code{\link{hclust}}}
}

