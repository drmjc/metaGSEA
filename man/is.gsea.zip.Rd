\name{is.gsea.zip}
\alias{is.gsea.zip}
\title{Is the argument a GSEA zip file?}
\usage{
  is.gsea.zip(x)
}
\arguments{
  \item{x}{a vector of file paths}
}
\value{
  a logical vector, same length as \code{x}
}
\description{
  Is the argument a GSEA zip file?
}
\examples{
f <- "my.zip"
is.gsea.zip(f)
is.gsea.zip(c(f,f,f,"/home"))
}
\author{
  Mark Cowley, 2011-08-24
}

