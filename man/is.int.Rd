\name{is.int}
\alias{is.int}
\title{Checks which elements of a vector are integers.}
\usage{
  is.int(x)
}
\arguments{
  \item{x}{a vector. usually this will be a numeric vector,
  but can also be a character}
}
\value{
  a logical vector, same length as x, with TRUE if each
  element is an integer
}
\description{
  This differs from is.integer which only checks the class
  of the entire vector.
}
\author{
  Mark Cowley, 2009-11-30
}

