\name{list2df}
\alias{list2df}
\title{Convert list of vectors to data.frame
Take a ragged \code{list} (ie a list of \code{vector}s of different lengths) and return a
\code{data.frame}, filling the shorter vectors with \code{NA}'s}
\usage{
  list2df(x)
}
\arguments{
  \item{x}{a list of vectors}
}
\value{
  a \code{data.frame}, filled with \code{NA}'s where
  necessary
}
\description{
  Convert list of vectors to data.frame Take a ragged
  \code{list} (ie a list of \code{vector}s of different
  lengths) and return a \code{data.frame}, filling the
  shorter vectors with \code{NA}'s
}
\author{
  Mark Cowley, 25/2/08
}

