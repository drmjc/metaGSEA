\name{gsea2adjmat}
\alias{gsea2adjmat}
\title{Convert a GSEA object into an adjacency matrix at the gene-level.}
\usage{
  gsea2adjmat(x)
}
\arguments{
  \item{x}{a GSEA list}
}
\value{
  a data.frame of 0 or 1 depending on whether each gene was
  found in each element of the GSEA list.
}
\description{
  From a GSEA object, extract the genes from the leading
  edges of the genesets, and then convert into an adjacency
  matrix.
}
\author{
  Mark Cowley, 2009-10-29
}
\seealso{
  \code{\link[mjcbase]{list2adjmat}}
  \code{\link{gsealist2adjmat}}
}

