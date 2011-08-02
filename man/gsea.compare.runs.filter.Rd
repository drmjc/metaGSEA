\name{gsea.compare.runs.filter}
\alias{gsea.compare.runs.filter}
\title{Function to filter a set of GSEA run comparisons, and/or to re-order the...}
\usage{gsea.compare.runs.filter(x, fdr, mode=c("any", "all")[1], sort.by=c("none", "avgabs", "range",
    "deltaAB", 1))
}
\description{Function to filter a set of GSEA run comparisons, and/or to re-order the
rows.}
\details{This also provides extensive ability to re-sort the filtered rows, via \code{sort.by} 
and the \code{\link{gsea.compare.runs.sort}}?

If you want to just sort the rows, see \code{\link{gsea.compare.runs.sort}}, or set \code{fdr=NULL}.}
\value{x with only the rows that satisfy the filter.}
\author{Mark Cowley, 2009-03-23}
\arguments{\item{x}{a data.frame from \code{\link{import.gsea.compare.runs}}}
\item{fdr}{a single numeric value, such as 0.05 or NULL}
\item{mode}{should \dQuote{any}, or \dQuote{all} of the values satisfy the filter?}
\item{sort.by}{one of \dQuote{none}, \dQuote{avgabs}, \dQuote{range}, \dQuote{deltaAB}, or an integer. see \code{\link{gsea.compare.runs.sort}}}
}

