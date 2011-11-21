\name{.qplot.data}
\alias{.qplot.data}
\title{.qplot.data
This code is taken from qvalue::qplot (from R2.9.0) in order to export the required data for my
LimmaGP GenePattern module}
\usage{
  .qplot.data(qobj, rng = c(0, 0.1), smooth.df = 3,
    smooth.log.pi0 = FALSE, ...)
}
\arguments{
  \item{qobj}{a qvalue object.}

  \item{rng}{see \code{\link[qvalue]{qplot}}}

  \item{smooth.df}{see \code{\link[qvalue]{qplot}}}

  \item{smooth.log.pi0}{see \code{\link[qvalue]{qplot}}}
}
\value{
  a named list of 4 vectors: \item{\code{lambda}}{the
  lambda values (usually 0...0.9)} \item{\code{pi0}}{the
  estimated pi0 for each lambda} \item{\code{spline}}{the
  cubed spline values for each lambda}
  \item{\code{pi0.estimate}}{the final pi0 estimate (from
  the qobj itself). length 1}

  for example:\verbatim{ List of 4\cr $ lambda : num [1:20]
  0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 ...\cr $ pi0 :
  num [1:20] 1 0.975 0.961 0.942 0.92 ...\cr $ spline : num
  [1:20] 0.987 0.97 0.953 0.936 0.92 ...\cr $ pi0.estimate:
  num 0.793\cr }
}
\description{
  .qplot.data This code is taken from qvalue::qplot (from
  R2.9.0) in order to export the required data for my
  LimmaGP GenePattern module
}
\author{
  Mark Cowley, 2009-11-30
}
\seealso{
  \code{\link[qvalue]{qplot}}
}

