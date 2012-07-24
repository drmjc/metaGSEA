\name{plot_gsea.leadingedge}
\alias{plot_gsea.leadingedge}
\title{Combination GSEA Leading Edge plot.}
\usage{
  plot_gsea.leadingedge(x, main)
}
\arguments{
  \item{x}{a GSEA list}

  \item{main}{the main plotting title}
}
\value{
  none. creates multiple plots.
}
\description{
  Combination GSEA Leading Edge plot.
}
\note{
  it's only worth using on filtered GSEA data - see
  \code{gsea.filter(*, N=50, direction="up")}
}
\author{
  Mark Cowley, 2009-10-29
}
\seealso{
  \code{\link{gsea.filter}}
  \code{\link{plot_gsea.leadingedge.HCL}}
  \code{\link{plot_gsea.leadingedge.heatmap}}
  \code{\link{plot_gsea.leadingedge.barplot}}
}

