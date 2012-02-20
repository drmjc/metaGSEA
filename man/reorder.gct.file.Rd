\name{reorder.gct.file}
\alias{reorder.gct.file}
\title{Reorder the rows of a GCT file.}
\usage{
  reorder.gct.file(gct.file, gct.file.out,
    method = c("var", "mean", "median", "max", "sum"),
    reverse = FALSE)
}
\arguments{
  \item{gct.file}{the path to a gct file}

  \item{gct.file.out}{the name of the GCT file to be
  created with the reordered rows}

  \item{method}{see details. Default = 'var'}

  \item{reverse}{if \code{TRUE}, then reverse the default
  sort order (see details for default sort order)}
}
\value{
  creates a file at gct.file.out, where the rows have been
  reordered
}
\description{
  see \code{\link{reorder.gct}} for more details
}
\author{
  Mark Cowley, 2011-03-16
}
\seealso{
  \code{\link{reorder.gct}}
}

