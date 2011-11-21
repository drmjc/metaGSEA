\name{import.gsea.cls}
\alias{import.gsea.cls}
\title{Import a cls file}
\usage{
  import.gsea.cls(file, as.factor = TRUE,
    enforce.zero.start = TRUE)
}
\arguments{
  \item{file}{the path to a valid cls file.}

  \item{as.factor}{logical: return a character, or factor?}

  \item{enforce.zero.start}{logical: check that the lowest
  class is 0. See details.}
}
\value{
  if cls file is continuous, either a character vector, or
  factor, one element per sample, of if continuous, a
  numeric is returned.
}
\description{
  Import a GenePattern cls file into R. CLS files can
  either be categorical, or continuous & both can be
  handled by this code.
}
\details{
  The cls file spec says for catgegorical cls files, values
  should start from 0. However some tools (eg SubMap)
  require them to start from 1. \code{enforce.zero.start}
  allows control over these situations.
}
\author{
  Mark Cowley, 2009-11-26
}
\seealso{
  \code{\link{import.gsea.gct}}
  \code{\link{import.gsea.clm}}
  \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#cls}
}

