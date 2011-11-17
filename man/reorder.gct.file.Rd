\name{reorder.gct.file}
\alias{reorder.gct.file}
\title{Reorder the rows of a GCT file.}
\usage{
  reorder.gct.file(gct.file, gct.file.out, method =
  c("var", "mean", "median", "sum"), reverse = FALSE)
}
\arguments{
  \item{gct.file}{see reorder.gct}

  \item{method}{see reorder.gct}

  \item{reverse}{see reorder.gct}

  \item{gct.file.out}{the name of the GCT file to be
  created with the reordered rows}
}
\value{
  creates a file at gct.file.out, where the rows have been
  reordered
}
\description{
  Details: see reorder.gct
}
\author{
  Mark Cowley, 2011-03-16
}

