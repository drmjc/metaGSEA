\name{plot.gsea.leadingedge.heatmap.auto}
\alias{plot.gsea.leadingedge.heatmap.auto}
\title{Generate a heatmap of genesets from a GSEA run, across a number of standard thresholds
The default settings will generate 4 heatmaps, of the up/down, top 50, and FDR<0.25 genesets.}
\usage{
  plot.gsea.leadingedge.heatmap.auto(x, file = NULL,
    N = c(50), FDR = c(0.25), P = NULL, FWER = NULL,
    pdf.width = 12, pdf.height = 12)
}
\arguments{
  \item{x}{a GSEA obect}

  \item{file}{the path to a PDF file, or \code{NULL} to use
  \code{dev.cur}}

  \item{N}{an integer vector of thresholds on the top N
  genesets, sorted by |NES|}

  \item{FDR}{an integer vector of FDR thresholds}

  \item{P}{an integer vector of unadjusted P-value
  thresholds}

  \item{FWER}{an integer vector of FWER thresholds}

  \item{pdf.width}{the width of the pdf device in inches}

  \item{pdf.height}{the height of the pdf device in inches}
}
\value{
  nothing. creates a pdf file.
}
\description{
  Generate a heatmap of genesets from a GSEA run, across a
  number of standard thresholds The default settings will
  generate 4 heatmaps, of the up/down, top 50, and FDR<0.25
  genesets.
}
\author{
  Mark Cowley, 2010-10-14
}

