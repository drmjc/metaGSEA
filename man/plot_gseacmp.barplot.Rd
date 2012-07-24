\name{plot_gseacmp.barplot}
\alias{plot_gseacmp.barplot}
\title{Taking a gseacmp object (see gsea.compare.runs.1gmt), create barplots of
the NES scores for each geneset}
\usage{
  plot_gseacmp.barplot(x, fdr.thresh = 0.05, do.par = TRUE,
    sub = "", max.label.length = 13, col = "lightblue",
    las = 2, nrow = 3, ncol = 4, legend.pos = "none",
    label.bars = c("stars", "FDR", "none"),
    label.fun = pvalue.stars)
}
\arguments{
  \item{x}{a gseacmp object (see gsea.compare.runs.1gmt)}

  \item{fdr.thresh}{Only keep genesets where at least 1
  result is < this threshold.}

  \item{do.par}{setup the plotting region? defaults to
  TRUE}

  \item{sub}{see ?title}

  \item{max.label.length}{max number of characters in a
  name. Middle characters will be replaced with ... to make
  up a name of at most max.label.length. Hint: choose an
  odd number.}

  \item{col}{the barplot's background colour. ?see barplot}

  \item{las}{see ?par}

  \item{nrow}{how many plots per page. Passed to par(mfrow)
  if do.par=TRUE}

  \item{ncol}{how many plots per page. Passed to par(mfrow)
  if do.par=TRUE}

  \item{legend.pos}{Where to put the legend ONLY on the
  first plot. Try "top" or "bottom"}

  \item{label.bars}{stars: statistically significant bars
  are labelled using stars. See label.fun}

  \item{FDR}{the FDR value is written on the bar, only if
  FDR < 0.25}

  \item{none}{no FDR labelling is done.}

  \item{label.fun}{Function used to work out how many stars
  to plot for each bar. Defaults to pvalue.stars.}
}
\description{
  Taking a gseacmp object (see gsea.compare.runs.1gmt),
  create barplots of the NES scores for each geneset
}
\author{
  Mark Cowley, 2009-03-23
}

