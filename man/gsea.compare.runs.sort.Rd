\name{gsea.compare.runs.sort}
\alias{gsea.compare.runs.sort}
\title{Function to sort the rows of a GSEA comparison object}
\usage{
  gsea.compare.runs.sort(x,
    sort.by = c("none", "avgabs", "range", "deltaAB", 1))
}
\arguments{
  \item{x}{a data.frame from
  \code{\link{import.gsea.compare.runs}}}

  \item{sort.by}{one of \dQuote{none}, \dQuote{avgabs},
  \dQuote{range}, \dQuote{deltaAB}, or an integer. see
  details}
}
\value{
  x with the rows resorted
}
\description{
  This provides extensive ability to re-sort the rows, via
  \code{sort.by}.\cr \dQuote{avgabs} means sort by
  decreasing avg absolute NES, so genesets with large NES's
  in all comparisons will be prioritised;\cr
  \dQuote{range}: genesets with large swings in NES across
  \emph{all} GSEA runs will be prioritised\cr
  \dQuote{deltaAB}: genesets are ranked by the magnitude of
  the difference between \emph{only} run 1 and 2; or\cr
  supply a single integer to specify which of the N runs
  you want to sort by.
}
\author{
  Mark Cowley, 2011-07-19
}
\seealso{
  \code{\link{gsea.compare.runs.filter}}
}

