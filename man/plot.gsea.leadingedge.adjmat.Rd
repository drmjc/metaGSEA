\name{plot.gsea.leadingedge.adjmat}
\alias{plot.gsea.leadingedge.adjmat}
\title{Generate a binary adjacency matrix image (like a heatmap with 2 colours).}
\usage{
  plot.gsea.leadingedge.adjmat(x,
    main = "Leading Edge Adjacency Matrix", cluster = TRUE)
}
\arguments{
  \item{x}{a GSEA object}

  \item{cluster}{if TRUE then do a 2-way hierarchical
  clustering, else order genes by decreasing frequency}

  \item{main}{the title}
}
\value{
  none. creates a plot.  Todo: introduce colour by looking
  up stat in the rnk object. handle the case of lots of
  genes nicer.
}
\description{
  The occurence of genes in each geneset is represented by
  blue squares. Like in the Broad's GSEA Leading Edge
  Viewer tool, genes are in columns, gensets are in rows.
}
\author{
  Mark Cowley, 2009-10-29
}

