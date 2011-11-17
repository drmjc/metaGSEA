\name{import.gsea.gct}
\alias{import.gsea.gct}
\title{Import a gct file}
\usage{
  import.gsea.gct(f, clm.file = NULL, check.names = FALSE)
}
\arguments{
  \item{f}{[character] the path to a gct file}

  \item{clm.file}{[optional] path to a clm.file}

  \item{check.names}{logical: R likes to make sure that the
  column names are 'syntactically valid', If you don't
  care, then choose FALSE. default is FALSE.}
}
\value{
  a \code{data.frame} of usually gene expression data.
  columns 1&2 are annotation, 3+ are data.
}
\description{
  Import a GenePattern gct file into R as a data.frame
  object. If a clm.file is also provided, then the columns
  of the gct are optionally subsetted, reordered, and
  renamed as per the clm.
}
\examples{
\dontrun{
f <- "ftp://ftp.broadinstitute.org/pub/genepattern/datasets/all_aml/all_aml_test.gct"
import.gsea.gct(f)
}
}
\author{
  Mark Cowley, 2009-07-27
}
\seealso{
  \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#gct},
  \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#clm}
}

