\name{reorder.gct}
\alias{reorder.gct}
\title{Function to reorder the rows of a GCT object}
\usage{
  reorder.gct(gct,
    method = c("var", "mean", "median", "max", "sum"),
    reverse = FALSE)
}
\arguments{
  \item{gct}{a gct object}

  \item{method}{see details. Default = 'var'}

  \item{reverse}{if \code{TRUE}, then reverse the default
  sort order (see details for default sort order)}
}
\description{
  take the data in a gct file & reorder the rows by
  different criteria:\cr - var: sort by variance, from high
  to low\cr - mean: sort by average abundance, from high to
  low\cr - median: sort by median abundance, from high to
  low\cr - max: sort by max abundance, from high to low\cr
  - sum: sort by sum, from high to low (useful if the GCT
  is an adjacency matrix/matrix of counts, or scores)\cr
}
\author{
  Mark Cowley, 2011-03-16
}
\seealso{
  \code{\link{reorder.gct.file}}
}

