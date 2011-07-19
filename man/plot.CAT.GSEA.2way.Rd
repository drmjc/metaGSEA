\name{plot.CAT.GSEA.2way}
\alias{plot.CAT.GSEA.2way}
\title{Plot a CAT plot upon 2 GSEA objects.}
\usage{plot.CAT.GSEA.2way(x, y, col=1, proportion=FALSE, add=FALSE, FDRthresh=1.01, main="")
}
\description{Plot a CAT plot upon 2 GSEA objects.}
\details{This is a helper function, you should probably consider using \code{\link{plot.CAT.GSEA}}.}
\value{none. makes a plot}
\author{Mark Cowley, 2011-07-19}
\arguments{\item{x}{a GSEA object}
\item{y}{a GSEA object}
\item{col}{the plot character foreground colour. default = 1 = \dQuote{black}}
\item{proportion}{if \code{TRUE}, x-axis is proportion of all genesets, if \code{FALSE}
(default), the number of genesets is plotted.}
\item{add}{logical: add this to existing plot? if \code{FALSE}, then create a new plotting device.}
\item{FDRthresh}{If you'd like the \code{pch} to change depending on which of A or
B is significant, then choose a threshold of 0.25 for example. Pch will
change depending on whether A and B are both significant, only one is
significant, or neither are.}
\item{main}{the plot title. set to \dQuote{} to ignore}
}

