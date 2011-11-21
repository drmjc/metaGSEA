\name{plclust_in_colour}
\alias{plclust_in_colour}
\title{modifiction of plclust for plotting hclust objects *in colour*!}
\usage{
  plclust_in_colour(hclust, lab = hclust$labels,
    lab.col = rep(1, length(hclust$labels)), hang = 0.1,
    ...)
}
\arguments{
  \item{hclust}{an hclust object}

  \item{lab}{a character vector of labels of the leaves of
  the tree}

  \item{lab.col}{colour for the labels; \code{NA}=default
  device foreground colour}

  \item{hang}{see \code{\link{plclust}}}

  \item{\dots}{args passed to \code{\link{text}}}
}
\value{
  none. Generates A display of hierarchical cluster with
  coloured leaf labels.
}
\description{
  modifiction of plclust for plotting hclust objects *in
  colour*!
}
\author{
  Eva Chan, 2009
}

