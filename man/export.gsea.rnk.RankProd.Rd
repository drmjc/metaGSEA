\name{export.gsea.rnk.RankProd}
\alias{export.gsea.rnk.RankProd}
\title{Take the results from a RankProd analysis, and produce a rnk file for GSEA.}
\usage{
  export.gsea.rnk.RankProd(rp, names = NULL, file,
    method = c("inv.log", "inverse"))
}
\arguments{
  \item{rp}{a rankprod result object.}

  \item{names}{a character vector of probe names, or if
  \code{NULL}, the rownames of x will be used, else error}

  \item{file}{the output rnk file name}

  \item{method}{one of \sQuote{inv.log}, or
  \sQuote{inverse}. see
  \code{\link{RankProd.signed.score}}}
}
\value{
  write out a rnk file, and invisibly return the rp score.
}
\description{
  The score is the inverse of the RP score, multiplied by
  100 or -100 for the up and down-regulated genes
  respectively. This typically produces scores in 3 to -3
  range
}
\author{
  Mark Cowley, 2009-01-09
}

