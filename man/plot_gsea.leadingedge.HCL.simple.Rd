\name{plot_gsea.leadingedge.HCL.simple}
\alias{plot_gsea.leadingedge.HCL.simple}
\title{Take a GSEA object and calculate an HCL between the leading edge genes.
Don't worry about the fancy geneset renaming that plot_gsea.leadingedge.HCL
gives you.}
\usage{
  plot_gsea.leadingedge.HCL.simple(x,
    main = "Leading edge similarities", xlab = "",
    hclust.method = "complete", ...)
}
\arguments{
  \item{x}{a GSEA object}

  \item{main}{the plot title}

  \item{xlab}{the x label}

  \item{hclust.method}{the hierarchical clustering method.
  see \code{\link{hclust}}}

  \item{\dots}{further arguments passed to
  \code{\link{plclust.gsea}}}
}
\value{
  none. creates an HCL dendrogram plot
}
\description{
  Take a GSEA object and calculate an HCL between the
  leading edge genes. Don't worry about the fancy geneset
  renaming that plot_gsea.leadingedge.HCL gives you.
}
\author{
  Mark Cowley, 2009-09-01
}
\seealso{
  \code{\link{plot_gsea.leadingedge.HCL}}
  \code{\link{plclust.gsea}}
  \code{\link{gsea.leadingedge.distance}}
}

