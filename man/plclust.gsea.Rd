\name{plclust.gsea}
\alias{plclust.gsea}
\title{Plot an hclust object, that has been made from GSEA genesets.}
\usage{
  plclust.gsea(hc, main = "Leading edge similarities",
    xlab = "", sub = NULL,
    yaxis = c("similarity", "distance"), lab.col = "black",
    hang = 0.03, cex = NULL, cex.main = 1, ...)
}
\arguments{
  \item{hc}{an object of class \code{hclust}}

  \item{main}{see \code{\link{par}}}

  \item{xlab}{see \code{\link{par}}}

  \item{sub}{see \code{\link{par}}}

  \item{yaxis}{what type of yaxis scale do you want?
  \dQuote{similarity}, or \dQuote{distance}?}

  \item{lab.col}{argument to change the colour of the text
  labels. (not yet implemented)}

  \item{hang}{see \code{\link{plclust}} default=0.03}

  \item{cex}{see \code{\link{par}}. default=NULL}

  \item{cex.main}{see \code{\link{par}}. default=1}

  \item{\dots}{arguments passed to \code{\link{plclust}}'}
}
\value{
  none. creates a plot
}
\description{
  Plot an hclust object, that has been made from GSEA
  genesets.
}
\author{
  Mark Cowley, 2009-09-04
}
\seealso{
  \code{\link{plclust}}
}

