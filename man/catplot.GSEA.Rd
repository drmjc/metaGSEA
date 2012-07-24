\name{catplot.GSEA}
\alias{catplot.GSEA}
\title{Function to generate CAT plots for >= 2 GSEA comparisons}
\usage{
  catplot.GSEA(..., gsea.list = NULL, labels = LETTERS,
    legend.pos = "bottomright", proportion = FALSE,
    FDRthresh = 1.01, main = "",
    method = c("1-vs-all", "all-pairs"))
}
\arguments{
  \item{\dots}{at least 2 GSEA objects, or:}

  \item{gsea.list}{an optional list of >= 2 GSEA objects.
  This overrides \dots}

  \item{labels}{the name of of each GSEA run. must be same
  length as \dots, or gsea.list.}

  \item{legend.pos}{where to add the names legend. default
  = \dQuote{bottomright}, set to \dQuote{none} to disable}

  \item{proportion}{if \code{TRUE}, x-axis is proportion of
  all genesets, if \code{FALSE} (default), the number of
  genesets is plotted.}

  \item{main}{the plot title.}

  \item{FDRthresh}{If you'd like the \code{pch} to change
  depending on which of A or B is significant, then choose
  a threshold of 0.25 for example. Pch will change
  depending on whether A and B are both significant, only
  one is significant, or neither are.}

  \item{method}{\dQuote{1-vs-all} means compare the first
  GSEA run to all subseqent runs; \dQuote{all-pairs} plots
  all pairwise comparisons of genesets.}
}
\description{
  This function expects all GSEA runs to be generated vs
  the \emph{same} gmt file, and with no filtering done.
  Aside from the filtering upon geneset size that GSEA
  performs, there should be a large overlap in terms of
  genesets
}
\author{
  Mark Cowley, 2008-12-08 (edited date)
}

