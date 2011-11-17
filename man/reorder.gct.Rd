\name{reorder.gct}
\alias{reorder.gct}
\title{Function to reorder the rows of a GCT file}
\usage{
  reorder.gct(gct, method = c("var", "mean", "median",
  "sum"), reverse = FALSE)
}
\arguments{
  \item{gct}{a gct object}

  \item{method}{see details. Defualt = 'var'}

  \item{reverse}{if TRUE, the reverse the default sort
  order (see details for default sort order)}
}
\description{
  take the data in a gct file & reorder the rows by
  different criteria: - var: sort by variance, from high to
  low - mean: sort by average abundance, from high to low -
  median: sort by median abundance, from high to low - sum:
  sort by sum, from high to low (useful if the GCT is an
  adjacency matrix/matrix of counts, or scores)
}
\author{
  Mark Cowley, 2011-03-16
}

