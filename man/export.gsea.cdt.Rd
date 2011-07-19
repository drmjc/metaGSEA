\name{export.gsea.cdt}
\alias{export.gsea.cdt}
\title{Cluster leading edge genesets, and export in Stanford cdt/gtr format.}
\usage{export.gsea.cdt(x, file, gtr=TRUE, rename=FALSE)
}
\description{Cluster leading edge genesets, and export in Stanford cdt/gtr format.}
\details{cdt files = Clustered Data Table; gtr files = gene tree file. These files
describe the ordering of usually a matrix of gene expression data, however
we have used this notation to represent the relationships between GSEA
leading edge genesets. You can view these files using the
\code{HierarchicalClusteringViewer} tool in \code{GenePattern}.

The order of the tree should be identical to \code{\link{plot.gsea.leadingedge.HCL}}, so
you could compare a static pdf image (from \code{\link{plot.gsea.leadingedge.HCL}}), or
interrogate a dynamically updating image using the \code{GenePattern} tool.}
\author{Mark Cowley, 2009-10-12}
\seealso{\code{\link{plot.gsea.leadingedge.HCL}}}
\references{Stanford File formats:
\url{http://smd.stanford.edu/help/formats.shtml} 
GenePattern file formats:
\url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html}}
\keyword{manip}
\keyword{IO}
\keyword{file}
\value{nothing. Creates a cdt file, and optionally a gtr file.}
\arguments{\item{x}{a gsea object (not the leading.edge)}
\item{file}{the path to a cdt file. This should have the .cdt extension}
\item{gtr}{if TRUE, then a gtr file will also be made. It gets its filename
from the cdt. if \code{FALSE}, then it will not be made (can't see why that
would be useful to you though?)}
\item{\dots}{Unprocessed arguments for export.gsea.gtr.}
\item{rename}{logical: rename the genesets to at most 40 characters? see \code{\link{gsea.rename.genesets}}}
}

