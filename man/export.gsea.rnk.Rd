\name{export.gsea.rnk}
\alias{export.gsea.rnk}
\title{Export a vector of values as a GSEA pre-ranked list.}
\usage{export.gsea.rnk(values, names, file, resort=TRUE, sigfigs=4)
}
\description{Export a vector of values as a GSEA pre-ranked list.}
\value{<none>, makes a file}
\author{Mark Cowley, 2008-09-08}
\arguments{\item{values}{a numeric vector}
\item{names}{the names associated with each value. If \code{NULL}, then the
\code{names(values)} are used}
\item{file}{the output file name. No hyphens, or a warning will be issued.}
\item{resort}{logical: re-sort the values such that they are decreasing? see order(,
decreasing=TRUE)}
\item{sigfigs}{the number of decimal places, defaults to 4.}
}

