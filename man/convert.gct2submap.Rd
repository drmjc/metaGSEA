\name{convert.gct2submap}
\alias{convert.gct2submap}
\title{Prepare data for use in SubMap.}
\usage{
  convert.gct2submap(gct.file, cls.file, chip.file,
    gct.out, cls.out,
    collapse.mode = c("max.avg", "mean", "median", "max.cv"),
    cv.thresh = NULL, topN.thresh = NULL)
}
\arguments{
  \item{gct.file}{standard input files to be modified.}

  \item{cls.file}{standard input files to be modified.}

  \item{chip.file}{used to map probe ID's to gene symbols}

  \item{gct.out}{output files.}

  \item{cls.out}{output files.}

  \item{collapse.mode}{how to choose the value for each
  gene.}

  \item{max.avg}{calculate the average for each probeset &
  choose the values from the probes with largest average
  value}

  \item{mean}{for N probesets mapping to the same gene,
  calculate the mean or median of those N values for each
  sample.}

  \item{median}{for N probesets mapping to the same gene,
  calculate the mean or median of those N values for each
  sample.}

  \item{max.cv}{calculate the CV for each probe & choose
  the values from the probe with highest cv}

  \item{cv.thresh}{[optional] CV threshold to use. CV's are
  percentages, so 50, or 100 is a good place to start.
  Thresholding is done after collapsing to unique genes.
  This overrides topN.thresh.}

  \item{topN.thresh}{[optional] take the top N genes when
  ranked by their CV's. Must leave cv.thresh=NULL.
  Thresholding is done after collapsing to unique genes.s}
}
\value{
  creates a new gct & cls file. The gct file has 1 row per
  gene which passed the filters & has gene symbols in
  column 1 (the Name column). The cls file starts from '1'
  instead of '0' which means it may not work with other
  GenePattern modules!
}
\description{
  Converts a gct file and cls file so that: - one row per
  gene symbol & gene symbol becomes the Name of each row,
  rather than the probe id -- different methods of
  collapsing rows - cls file starts from 1, rather than 0 -
  choose the most variable genes, either by a coefficient
  of variation threshold (eg 50) or the top N most variable
  genes, as defined by CV.
}
\author{
  Mark Cowley, 2010-07-06
}

