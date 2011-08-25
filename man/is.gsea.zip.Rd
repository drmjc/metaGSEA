\name{is.gsea.zip}
\alias{is.gsea.zip}
\title{Is the argument a GSEA zip file?}
\usage{is.gsea.zip(x)
}
\description{Is the argument a GSEA zip file?}
\value{a logical vector, same length as \code{x}}
\author{Mark Cowley, 2011-08-24}
\arguments{\item{x}{a vector of file paths}
}
\examples{f <- "my.zip"
is.gsea.zip(f)
is.gsea.zip(c(f,f,f,"/home"))}

