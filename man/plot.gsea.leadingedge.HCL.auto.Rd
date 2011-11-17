\name{plot.gsea.leadingedge.HCL.auto}
\alias{plot.gsea.leadingedge.HCL.auto}
\title{Automated method for doing a hierarchical clustering of genesets, based on
each sets leading edge genes, for all gsea results in a list of gsea
results, over a range of different best N, FDR, P, FWER thresholds, and for
up/down genesets.}
\usage{
  plot.gsea.leadingedge.HCL.auto(x, file = NULL, N = c(50),
  FDR = c(0.25), P = NULL, FWER = NULL, max.distance = 1)
}
\arguments{
  \item{x}{either a gsea list, or a list of gsea lists. see
  import.gsea}

  \item{file}{either the path to a pdf file which will be
  created, or NULL to use the current plotting device}

  \item{N}{vectors of 0 or more thresholds to use. eg
  N=c(50, 100) means plot the top 50, or the}

  \item{FDR}{vectors of 0 or more thresholds to use. eg
  N=c(50, 100) means plot the top 50, or the}

  \item{P}{vectors of 0 or more thresholds to use. eg
  N=c(50, 100) means plot the top 50, or the}

  \item{FWER}{vectors of 0 or more thresholds to use. eg
  N=c(50, 100) means plot the top 50, or the top 100
  genesets}

  \item{max.distance}{(experimental) threshold to filter
  out genesets that are poorly connected to each other.
  This is a distance threshold, so 1.0 means all genesets,
  0.9 is a good place to start.}
}
\value{
  makes a number of plots to a pdf file, or the current
  plotting device See also: plot.gsea.leadingedge.HCL,
  gsea.leadingedge.distance
}
\description{
  Automated method for doing a hierarchical clustering of
  genesets, based on each sets leading edge genes, for all
  gsea results in a list of gsea results, over a range of
  different best N, FDR, P, FWER thresholds, and for
  up/down genesets.
}
\author{
  Mark Cowley, 2009-09-03
}

