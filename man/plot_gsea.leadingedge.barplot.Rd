\name{plot_gsea.leadingedge.barplot}
\alias{plot_gsea.leadingedge.barplot}
\title{GSEA Leading Edge analysis: Barplot of number of genesets each gene occured
in.}
\usage{
  plot_gsea.leadingedge.barplot(x, min.count = 0,
    xlab = "", ylab = "Number Of Gene Sets",
    main = "Leading Edge Barplot", horiz = FALSE,
    col = "blue", cex.names = 0.6, ...)
}
\arguments{
  \item{x}{a GSEA object. Usually this has been subset to a
  managable number of genesets}

  \item{min.count}{filter out genes found in only a few
  genesets. min.count=2 is useful in practice.}

  \item{xlab}{see \code{\link{plot}}}

  \item{ylab}{see \code{\link{plot}}}

  \item{main}{see \code{\link{plot}}}

  \item{horiz}{logical: if TRUE, then the bars will be
  horizontal; FALSE they are vertical (default).}

  \item{col}{the colour of the bars}

  \item{cex.names}{the character expansion factor of the
  gene names. see \code{\link{plot}}}

  \item{\dots}{additional arguments passed to barplt}
}
\value{
  invisibly returns the named vector of counts.
}
\description{
  This tool emulates the barplot found in the GSEA Leading
  Edge Tool. Typically, the user selects a subset of
  genesets of interest (use \code{\link{gsea.filter}}, or
  \code{\link{gsea.split}} to do this) and then wants to
  find genes that are present in lots of those genesets.
  This tool gives you a barplot indicating the number of
  genesets that each gene was found in.
}
\author{
  Mark Cowley, 2009-10-28
}
\seealso{
  \code{\link{barplot}}
}

