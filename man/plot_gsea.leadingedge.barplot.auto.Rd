\name{plot_gsea.leadingedge.barplot.auto}
\alias{plot_gsea.leadingedge.barplot.auto}
\title{Function to automatically plot barplot's of filtered GSEA data.}
\usage{
  plot_gsea.leadingedge.barplot.auto(x, file = NULL,
    N = c(50), FDR = c(0.25), P = NULL, FWER = NULL,
    min.count = 2, horiz = FALSE, col = "blue",
    cex.names = 0.6, write.delim = FALSE, ...)
}
\arguments{
  \item{x}{a GSEA list}

  \item{file}{the path to an output PDF file. if
  \code{NULL}, then plots are printed to current device.}

  \item{N}{if not \code{NULL}, then choose the top N most
  significant genesets}

  \item{FDR}{if not \code{NULL}, then choose those genesets
  with FDR < the threshold}

  \item{P}{if not \code{NULL}, then choose those genesets
  with nominal P < the threshold}

  \item{FWER}{if not \code{NULL}, then choose those
  genesets with FWER < the threshold}

  \item{min.count}{filter out genes found in only a few
  genesets. min.count=2 is useful in practice.}

  \item{horiz}{logical: horizontal or vertical plot?}

  \item{col}{the colour of the bars}

  \item{cex.names}{the character expansion factor of the
  gene names. see \code{\link{plot}}}

  \item{\dots}{arguments passed to
  \code{\link{plot_gsea.leadingedge.barplot}}}

  \item{write.delim}{logical: if \code{TRUE}, then write
  the data plotted as a tab delimited txt file, named
  \code{file + ".txt"}. Default is \code{FALSE}. If
  \code{file=NULL}, then write.delim will be set to
  \code{FALSE}.}
}
\value{
  none. Generates a number of plots, either to a pdf if
  \code{file!=NULL}, or \code{\link{dev.cur}}
}
\description{
  Function to automatically plot barplot's of filtered GSEA
  data.
}
\note{
  Usually, I filter by up or down, by top 50, and by
  FDR<0.25. This generates 4 plots
}
\author{
  Mark Cowley, 2010-10-14
}

