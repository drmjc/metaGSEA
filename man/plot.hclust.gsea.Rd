\name{plot.hclust.gsea}
\alias{plot.hclust.gsea}
\title{Plot an hclust object, that has been made from GSEA genesets.}
\usage{plot.hclust.gsea(hc, main="Leading edge similarities", xlab="", sub, yaxis=c("similarity", "distance"), lab.col="black", hang=0.03, cex, cex.main=1, ...)
}
\description{Plot an hclust object, that has been made from GSEA genesets.
This is a convenience wrapper around \code{\link{plclust.gsea}}}
\value{none. creates a plot}
\seealso{\code{\link{plclust}}}
\author{Mark Cowley, 2009-09-04}
\arguments{\item{hc}{an object of class \code{hclust}}
\item{main}{see \code{\link{par}}}
\item{xlab}{see \code{\link{par}}}
\item{sub}{see \code{\link{par}}}
\item{yaxis}{what type of yaxis scale do you want? \dQuote{similarity}, or \dQuote{distance}?}
\item{lab.col}{argument to change the colour of the text labels. (not yet
implemented)}
\item{hang}{see \code{\link{plclust}} default=0.03}
\item{cex}{see \code{\link{par}}. default=NULL}
\item{cex.main}{see \code{\link{par}}. default=1}
\item{\dots}{arguments passed to \code{\link{plclust}}'}
}

