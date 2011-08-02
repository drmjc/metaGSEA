\name{plot.gsea.leadingedge.heatmap}
\alias{plot.gsea.leadingedge.heatmap}
\title{Generate a heatmap of genesets from a GSEA run...}
\usage{plot.gsea.leadingedge.heatmap(x, main="Geneset similarity", gridcolour="lightgrey")
}
\description{Generate a heatmap of genesets from a GSEA run}
\details{Inspired by the heatmap plot in the BROAD's GSEA Leading Edge Viewer tool,
this function creates a plot
which represents the similarity/dis-similarity between genesets with
colours. The colour is based on the
\code{Jacquard}. see \code{\link{dist.gsea}}}
\note{Recommend saving to file using pdf(..., width=12, height=12)}
\value{nothing. generates a plot.}
\author{Mark Cowley, 2009-10-29}
\arguments{\item{x}{a GSEA object}
\item{main}{the plot title}
\item{gridcolour}{the colour of the grid to overlay on top of the image. 
see sepcolor parameter of \code{\link[gplots]{heatmap.2}}}
}

