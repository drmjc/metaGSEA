\name{export.gsea.gmt}
\alias{export.gsea.gmt}
\title{Export a gmt file...}
\usage{export.gsea.gmt(x, f, description=rep("na", length(x)))
}
\description{Export a gmt file
Genesets for BROAD's GSEA are usually represented either as \code{gmt} files, or
\code{gmx} files. These can be opened in \code{excel}, and the \code{gmt} contains the first row
being the geneset name, and subsequent rows containing genes within each
geneset. These genes do not have to be the same length in each geneset.}
\details{See the detailed file specification, and the \code{\link{import.gsea.gmt}}
function.}
\value{none. it creates a gmt file.}
\author{Mark Cowley, 2009-04-29}
\seealso{\code{\link{import.gsea.gmt}}}
\references{\url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#gmt}}
\keyword{file}
\keyword{IO}
\arguments{\item{x}{A list of genesets. \code{x} must have names which are the geneset names,
and the elements are character vectors.}
\item{f}{the file name}
\item{description}{a character vector for the description column. it is recycled if its too short.
The order must match up with the rows of \code{x}.}
}

