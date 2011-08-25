\name{list2df}
\alias{list2df}
\title{Convert list of vectors to data...}
\usage{list2df(x)
}
\description{Convert list of vectors to data.frame
Take a ragged \code{list} (ie a list of \code{vector}s of different lengths) and return a
\code{data.frame}, filling the shorter vectors with \code{NA}'s}
\value{a \code{data.frame}, filled with \code{NA}'s where necessary}
\author{Mark Cowley, 25/2/08}
\arguments{\item{x}{a list of vectors}
}

