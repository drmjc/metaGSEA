\name{export.broad.cls}
\alias{export.broad.cls}
\title{Export a cls file for GenePattern.}
\usage{export.broad.cls(labels, file, continuous=FALSE, start.value=0)
}
\description{Export a cls file for GenePattern.
Export a cls file for GenePattern. cls files represent the class labels - a
way of telling GenePattern which experimental group each array belongs to.
They can be discrete, or continuous.}
\details{In the discrete case, each cls can contain only only parameter, thus the
argument must be a vector of labels. If you'd like control over which is
gets assigned to 0, 1, 2... then use pass in a factor of your making.

In the continuous case, you can have >= 1 parameter.  either pass in a
numeric vector; or a matrix/data.frame where each row represents a
different parameter; or a list where each element is a vector, and each
element represents a parameter.}
\value{none. creates a cls file.}
\note{update 2009-07-08: cls labels can't have spaces in them.}
\author{Mark Cowley, 2009-06-26}
\seealso{\code{\link{import.gsea.cls}}}
\references{\url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#cls}}
\keyword{file}
\keyword{IO}
\arguments{\item{labels}{the labels. see description.}
\item{file}{the output file name.}
\item{continuous}{are the labels continuous or discrete?
\code{TRUE}/\code{FALSE}, respectively.}
\item{start.value}{most cls files start from 0, however you can change this if you like.}
}
\examples{cls <- c(rep("wt",5), rep("mut", 5))
f <- tempfile()
export.broad.cls(cls, f, FALSE, 0)
cls2 <- rnorm(10)
f2 <- tempfile()
export.broad.cls(cls2, f2, TRUE, 0)}

