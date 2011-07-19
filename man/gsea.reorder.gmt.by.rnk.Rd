\name{gsea.reorder.gmt.by.rnk}
\alias{gsea.reorder.gmt.by.rnk}
\title{Reorder the genes within a gmt object to match the order of genes within a...}
\usage{gsea.reorder.gmt.by.rnk(gmt, rnk)
}
\description{Reorder the genes within a gmt object to match the order of genes within a
rnk object}
\value{a gmt just like input, but with the genes reordered to match the
order inside the rnk}
\author{Mark Cowley, 2009-12-07}
\arguments{\item{gmt}{a gmt object (ie a named list of genes)}
\item{rnk}{a rnk object (ie named list of gene scores, where genes are the
names)}
}

