\name{export.gsea.odf}
\alias{export.gsea.odf}
\title{Export a table of data to a GenePattern odf file}
\usage{
  export.gsea.odf(data, file, header.lines = NULL,
  colclasses = NULL)
}
\arguments{
  \item{data}{a \code{\link{data.frame}} or
  \code{\link{matrix}}}

  \item{file}{the path to the output odf file}

  \item{header.lines}{a character vector of header lines,
  or \code{NULL}}

  \item{colclasses}{a character vector of column class
  names. if \code{NULL}, then they are aut-determined from
  the \code{colclass} of each column, and mapped to
  GenePattern column class names}
}
\value{
  none. writes a file
}
\description{
  ODF files are tab delimited documents with fairly
  extensive header lines. Strongly recommend you look at
  the documentation for odf files (see below), and to see
  it in action, see \code{\link{export.gsea.odf.lmFit}}.
}
\note{
  header lines.\cr the header lines should not contain "ODF
  1.0", "HeaderLines=X", though the code ignore them. The
  only required field that you should consider providing is
  the \dQuote{Model=DataSet}. The \dQuote{DataLines=XYZ}
  header line is determined automatically.
}
\author{
  Mark Cowley, 2009-11-30
}
\seealso{
  \code{\link{export.gsea.odf.lmFit}}
  \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#odf}
}

