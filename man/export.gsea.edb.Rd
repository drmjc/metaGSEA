\name{export.gsea.edb}
\alias{export.gsea.edb}
\title{Export a GSEA edb object to an \code{XML} file.}
\usage{
  export.gsea.edb(edb, file)
}
\arguments{
  \item{edb}{an edb object.}

  \item{file}{the name of an xml file. Usually this is
  inside a dir called edb, and ends in \dQuote{.edb}}
}
\value{
  creates an XML file
}
\description{
  Export a GSEA edb object to an \code{XML} file.
}
\author{
  Mark Cowley, 2009-10-12
}
\seealso{
  \code{\link{import.gsea.edb}}, \code{\link[XML]{saveXML}}
}
\keyword{file}
\keyword{IO}

