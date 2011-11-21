\name{export.gsea.chip}
\alias{export.gsea.chip}
\title{Export a GSEA chip file
Export a GSEA chip file -- this basically ensures that the correct column
names are exported.}
\usage{
  export.gsea.chip(x, f, resort.rows = TRUE,
    collapse.rows = TRUE, na.strings = "---")
}
\arguments{
  \item{x}{a matrix-like object. see details, and
  \code{\link{import.gsea.chip}}}

  \item{f}{the output file name}

  \item{resort.rows}{re-order the rows using \code{order}
  on column 1.}

  \item{collapse.rows}{ensure that each probe is only
  outputted once}

  \item{na.strings}{what character value to replace all of
  the NA's with. eg \dQuote{} or \dQuote{NULL} or
  \dQuote{---}}
}
\value{
  none - creates the file.
}
\description{
  Export a GSEA chip file -- this basically ensures that
  the correct column names are exported.
}
\details{
  Assumes that the user provided these columns in order: 1.
  Probe ID 2. Gene Symbol 3. Gene Description 4. optional
  Gene Symbol Aliases
}
\author{
  Mark Cowley, 2008-09-01
}
\references{
  \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#chip}
}
\seealso{
  \code{\link{import.gsea.chip}}
}
\keyword{file}
\keyword{IO}

