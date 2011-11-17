\name{gsea.compare.runs.1gmt}
\alias{gsea.compare.runs.1gmt}
\title{Compare N GSEA runs vs a single gmt file (eg lots of rnk files vs c2_all)}
\usage{
  gsea.compare.runs.1gmt(dirs, gsealist, outfile, method =
  c("intersect", "union")[1])
}
\arguments{
  \item{dirs}{a vector of gsea directories}

  \item{gsealist}{optional list of gsea runs that you've
  already imported. over-rides <dirs>}

  \item{outfile}{if specified, an XLS output file will be
  created.}

  \item{method}{which genesets to keep if they are not the
  same? \dQuote{intersect}: just the genesets found in all
  GSEA runs; \dQuote{union}: all the genesets found in any
  GSEA run. \code{NA}'s will be used for those genesets
  that were not observed.}
}
\value{
  a data.frame of \dQuote{geneset}, \dQuote{size},
  \dQuote{NES}, \dQuote{FDR}, \dQuote{NES}, \dQuote{FDR},
  \dots
}
\description{
  Compare N GSEA runs vs a single gmt file (eg lots of rnk
  files vs c2_all)
}
\author{
  Mark Cowley, 2009-12-16
}

