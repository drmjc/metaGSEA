\name{export.broad.cls}
\alias{export.broad.cls}
\title{Export a cls file for GenePattern.
Export a cls file for GenePattern. cls files represent the class labels - a
way of telling GenePattern which experimental group each array belongs to.
They can be discrete, or continuous.}
\usage{
  export.broad.cls(labels, file, continuous = FALSE,
    start.value = 0, check.levels = TRUE)
}
\arguments{
  \item{labels}{the labels. see description.}

  \item{file}{the output file name.}

  \item{continuous}{are the labels continuous or discrete?
  \code{TRUE}/\code{FALSE}, respectively.}

  \item{start.value}{most cls files start from 0, however
  you can change this if you like.}

  \item{check.levels}{logical: if \code{TRUE}, then the
  order of the levels must match the order of the samples
  (ie avoid the ambiguous cls files )}
}
\value{
  none. creates a cls file.
}
\description{
  In the discrete case, each cls can contain only one
  parameter, thus the argument must be a vector of labels.
  If you'd like control over which class gets assigned to
  0, 1, 2... then use pass in a factor of your choosing &
  the levels of the factor will guide the class names.
}
\details{
  There is some ambiguity in the way samples are assigned
  to classes. eg:\cr 2 2 1 \cr # classA classB \cr 1 0 \cr
  The majority of tools would say the first sample belongs
  to class B (eg ScorebyClassComp), others (eg GSEA) would
  say that the first sample is classA, and maps all
  instances of a '1' to classA, and all instances of a '0'
  to the 2nd class. And others like PGSEA enforce the
  levels c(1,2). If \code{check.levels=TRUE}, and labels is
  a \code{factor}, then this checks that the order of the
  levels matches the order of the samples. In general, I
  try write non-ambiguous cls files, but sometimes it is
  unavoidable, in which case set \code{check.levels=FALSE}.

  In the continuous case, you can have >= 1 parameter.
  either pass in a numeric vector; or a matrix/data.frame
  where each row represents a different parameter; or a
  list where each element is a vector, and each element
  represents a parameter.
}
\note{
  update 2009-07-08: cls labels can't have spaces in them.
}
\examples{
cls <- c(rep("wt",5), rep("mut", 5))
f <- tempfile()
export.broad.cls(cls, f, FALSE, 0)
cls2 <- rnorm(10)
f2 <- tempfile()
export.broad.cls(cls2, f2, TRUE, 0)

cls3 <- factor(c("ClassB", "ClassA"), levels=c("ClassA", "ClassB"))
export.broad.cls(cls3, f2, FALSE, 0, check.levels=TRUE)
export.broad.cls(cls3, f2, FALSE, 0, check.levels=FALSE)
}
\author{
  Mark Cowley, 2009-06-26
}
\references{
  \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#cls}
}
\seealso{
  \code{\link{import.gsea.cls}}
}
\keyword{file}
\keyword{IO}

