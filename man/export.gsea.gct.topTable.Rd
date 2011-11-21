\name{export.gsea.gct.topTable}
\alias{export.gsea.gct.topTable}
\title{Export a GCT file, from a GSEA topTable plus data table.
This is a companion function to \code{\link{export.gsea.rnk.topTable}} which makes a gct
file which can then be used for annotating GSEA output with heatmaps.}
\usage{
  export.gsea.gct.topTable(tt, rma, chip, file,
    missing = "", verbose = TRUE, ...)
}
\arguments{
  \item{tt}{a GSEA top-table at the probe-level}

  \item{rma}{a data.frame of expression data with rownames
  being the probesetID}

  \item{chip}{the chip file that will be used in the
  gseaPreRanked analysis. Either the path to the file name,
  or the imported version as a data.frame. see
  \code{\link{import.gsea.chip}}}

  \item{file}{the fully specified path for the gct file}

  \item{missing}{how would you like}

  \item{verbose}{logical: print messages?}

  \item{\dots}{further arguments passed to
  \code{\link{export.gsea.gct}}}
}
\description{
  Export a GCT file, from a GSEA topTable plus data table.
  This is a companion function to
  \code{\link{export.gsea.rnk.topTable}} which makes a gct
  file which can then be used for annotating GSEA output
  with heatmaps.
}
\author{
  Mark Cowley, 2009-08-13
}
\seealso{
  \code{\link{export.gsea.gct}}
}

