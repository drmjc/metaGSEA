\name{plot.gsea.leadingedge.HCL.1clusterPerPage}
\alias{plot.gsea.leadingedge.HCL.1clusterPerPage}
\title{Take a GSEA object, do an HCL between the leading edge genes, cut the HCL,...}
\usage{plot.gsea.leadingedge.HCL.1clusterPerPage(x, main="Leading edge similarities", xlab="", hclust.method="complete",
    h=0.99, ...)
}
\description{Take a GSEA object, do an HCL between the leading edge genes, cut the HCL, 
then plot 1 HCL per page of a PDF.}
\author{Mark Cowley, 2009-09-01}
\arguments{\item{x}{a GSEA object}
\item{main}{the plot title}
\item{xlab}{the x label}
\item{hclust.method}{the hierarchical clustering method. see \code{\link{hclust}}}
\item{h}{numeric scalar or vector with heights where the tree should be
cut. (see cutree)}
\item{\dots}{see \code{\link{plot.gsea.leadingedge.HCL.simple}}}
}

