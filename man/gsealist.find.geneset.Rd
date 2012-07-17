\name{gsealist.find.geneset}
\alias{gsealist.find.geneset}
\title{From a list of GSEA's, find the elements that contain particular genesets
of interest.}
\usage{
  gsealist.find.geneset(gsea.list, patterns)
}
\arguments{
  \item{gsea.list}{a list of GSEA objects.}

  \item{patterns}{a vector of patterns which should be
  unique enough to find just one cluster}
}
\value{
  a vector of geneset indices where each pattern was found.
  IF one of the patterns is not unique enough, then a list
  is returned.
}
\description{
  From a list of GSEA's, find the elements that contain
  particular genesets of interest.
}
\section{Scenario}{
  You've used \code{\link{gsea.split}} to split a LARGE
  GSEA run into ~100 smaller, related GSEA runs. From a
  small list of genesets of interest (represented by
  patterns), you want to know which of these smaller GSEA
  runs contain your genesets of interest.
}
\author{
  Mark Cowley, 2009-10-16
}

