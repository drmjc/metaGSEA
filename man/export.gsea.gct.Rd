\% Generated by roxygen2 (4.0.2): do not edit by hand
\name{export.gsea.gct}
\alias{export.gsea.gct}
\title{Function to export a GenePattern gct file.}
\usage{
export.gsea.gct(data, description = NULL, file = NULL, chip = NULL,
  round = 4, version = "#1.2", missing = "", ...)
}
\arguments{
\item{data}{a \code{matrix} or \code{data.frame} WITH rownames of all numeric
 values, or a \code{data.frame} from running \code{\link{import.gsea.gct}};
ie first 2 columns are \dQuote{Name}, \dQuote{Description}.}

\item{description}{a vector of annotations. IF it has names, then these
names must match the rownames of data, in which case we make sure they're
in the same order. If no names, then we assume that they're in the same
order. If \code{NULL}, then you must either provide a \code{chip}, or have a
\dQuote{Description} column in the data.}

\item{file}{the output file name}

\item{chip}{instead of specifying description, you can specify a GSEA chip
object, and a description will be made for you. This overrides
\code{description}. Default is \code{NULL} to ignore.}

\item{round}{the number of digits to round the numbers to - default=4}

\item{version}{The GCT file version, to go in the first line.}

\item{missing}{the string to use for missing data (in the expression data).}

\item{\dots}{Currently unused.}
}
\description{
\code{data} can be either:\cr
a \code{matrix} or \code{data.frame} of numeric values,
in which case you can provide an optional \code{description} or \code{chip}
object to populate the \dQuote{Description} column;
or a \code{data.frame} containing \sQuote{Name}
and \dQuote{Description} columns (see \code{\link{import.gsea.gct}})
}
\examples{
\dontrun{
my.gct <- import.gsea.gct("./my.gct", "./my.cls")
export.gsea.gct(data=my.gct, file="my2.gct")
}
}
\author{
Mark Cowley, 2008-08-07
}

