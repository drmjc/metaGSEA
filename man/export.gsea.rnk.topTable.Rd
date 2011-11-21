\name{export.gsea.rnk.topTable}
\alias{export.gsea.rnk.topTable}
\title{Export a limma topTable to a GSEA pre-ranked list.
The 'rnk' file contains a column of ID's, and a column of sorted values, and
is typically used as input to GSEA in pre-ranked mode.}
\usage{
  export.gsea.rnk.topTable(tt, file,
    values = c("t", "logFC", "P.Value"), probe2gene = NULL,
    convert2symbol = FALSE, verbose = TRUE)
}
\arguments{
  \item{tt}{a limma \code{\link[limma]{topTable}}}

  \item{file}{the output file <****.rnk>}

  \item{values}{choose the column name that contains the
  values. This must be in the colnames of the topTable. If
  P.Value is chosen, then the -log10 P will be calculated,
  thus small P-values become large +ve numbers. Default is
  \dQuote{t}.}

  \item{probe2gene}{a 2 column matrix-like object with
  mappings from probes 2 genes. Preferably, this will be
  the chip file used for the GSEA. This allows the
  identification of the best performing probe per gene,
  rather than GSEA's "Max_Probe" method which chooses the
  probe with the larges value which is unsuitable for
  down-regulated genes where one poorly performing probe
  could be ~0. Leave as NULL to ignore this.}

  \item{convert2symbol}{logical: if
  \code{probe2gene!=NULL}, then after collapsing to the
  best probe-per-gene, export the genesymbol (\code{TRUE}),
  or the probe ID (\code{FALSE})?}

  \item{verbose}{logical: verbose messages}
}
\description{
  GSEA's Max_probe algorithm\cr GSEA has 2 modes for
  collapsing probes to genes: \dQuote{Max_probe} and
  \dQuote{Median_probe}, both work well when the input data
  is a GCT file, consisting of expression levels. However,
  I feel pretty strongly that GSEApreRanked needs a 3rd
  option which chooses the best probe, based on the one
  that obtained the largest t-statistic in either
  direction, a \dQuote{Best_probe} if you will. If you
  provide this function a \code{probe2gene} map (from probe
  ids to gene symbols), then we will automatically export
  just the best probes, ie 1 row per gene, still using the
  probe ID as the identifier (so that we can still use the
  same chip file).
}
\details{
  \code{probe2gene}\cr The justification for
  \code{probe2gene} is described above. \code{probe2gene}
  should be a 2 column \code{matrix-like} object with
  mappings from probes 2 genes. Preferably, this will be
  from the same chip file that will be used for running
  GSEA.
}
\author{
  Mark Cowley, 8/2/08
}
\seealso{
  \code{\link[limma]{topTable}}
}

