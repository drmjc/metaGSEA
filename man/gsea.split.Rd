\name{gsea.split}
\alias{gsea.split}
\title{Split a GSEA object based on leading edge similarities.}
\usage{gsea.split(x, dist=0.99, min.size=1)
}
\description{Split a GSEA object based on leading edge similarities.
Take a GSEA object upon G genesets, and split it into k GSEA objects,
where all genesets within each object has similarity to all other genesets.}
\details{Internally, this function first runs \code{\link{gsea.leadingedge.distance}}, then
\code{\link{hclust}}, and finally \code{\link{cutree}}. The resulting cluster assignments
are then used to subet GSEA objects into smaller GSEA objects with only
specified genesets within them. You can then export them, and load each
subset in the LeadingEdge Viewer tool.}
\value{a named list of GSEA objects, one per cluster of genesets. The
names are cluster001, cluster002, ...}
\author{Mark Cowley, 2009-10-16}
\arguments{\item{x}{a GSEA object}
\item{dist}{the distance to cut at. 0.99 is a large distance (close to the
root of the tree).}
\item{min.size}{discard geneset clusters that have fewer than min.size
genesets in them.}
}

