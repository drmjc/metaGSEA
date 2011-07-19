\name{ucounts}
\alias{ucounts}
\title{How many times does each term in x appear?}
\usage{ucounts(x, return.unique=TRUE, issorted=FALSE, sort=TRUE)
}
\description{How many times does each term in x appear?}
\details{\code{\link{table}} does a similar job, but does not provide the option
of \code{return.unique=FALSE}}
\value{a vector of counts, named by the elements of x}
\author{Mark Cowley, 11 May 2006}
\arguments{\item{x}{a character vector}
\item{return.unique}{logical: if TRUE, the numeric counts of the unique elements in x
are returned (the names indicate the elements in x). NB
the order of the returned elements is the same order that
unique(x) returns.

if FALSE, then return the numeric counts of each element
in x. Useful if you then want to work out which
elements have count > 1 and quickly view those elements.
NB, the order remains the same as x.}
\item{issorted}{logical: to improve the speed of this code, the algorithm first sorts the data.
If it's already sorted, then save time by setting \code{issorted=FALSE}.}
\item{sort}{logical: whether to sort the result vector from highest counts to lowest?}
}
\examples{ucounts(c("a", "a", "b", "c", "c"), TRUE)
# a b c
# 2 1 2
ucounts(c("a", "a", "b", "c", "c"), FALSE)
# a a b c c
# 2 2 1 2 2
ucounts(c("a", "a", "b", "c", "c"), FALSE) > 1
#    a     a     b     c     c
# TRUE  TRUE FALSE  TRUE  TRUE}

