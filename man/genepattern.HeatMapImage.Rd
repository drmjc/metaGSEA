\name{genepattern.HeatMapImage}
\alias{genepattern.HeatMapImage}
\title{From a GCT file, create a HeatMap Image, using the jar from Broad's GenePattern.
The Broad's HeatMapViewer tool generates nice looking, customisable heatmaps,
IMHO nicer than anything that \R can give, except maybe \code{\link[gplots]{heatmap.2}}.
If you have a bunch of GCT files that need converting to heatmap images fast,
then this is the code for you.}
\usage{
  genepattern.HeatMapImage(in.file, out.file = NULL, method
  = c("local", "server"), gp.connection = NULL, format =
  "png", row.size = 16, column.size = 16, show.grid =
  "yes", grid.color = "0:0:0", show.row.descriptions =
  "yes", show.row.names = "yes", rows.to.highlight = "",
  row.highlight.color = "255:0:0", color.scheme = "row
  normalized", color.palette = "", use.color.gradient =
  "no", memory = "-Xmx512m", verbose = FALSE, overwrite =
  FALSE)
}
\arguments{
  \item{in.file}{a single gct/res or odf file.}

  \item{out.file}{if NULL, the outfile will be in the same
  dir as in.file, with the extenstion changed dependant on
  format}

  \item{method}{\dQuote{local} or \dQuote{server}}

  \item{gp.connection}{if \code{method="server"} then you
  should specify a connection to the server. see examples.}

  \item{format}{picture format, one of: \dQuote{png},
  \dQuote{tiff}, \dQuote{jpeg}, \dQuote{eps}, \dQuote{bmp}.
  \dQuote{png} is default & recommended. \dQuote{tiff} look
  good for publication.}

  \item{row.size}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{column.size}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{show.grid}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{grid.color}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{show.row.descriptions}{see HeatMapImage help at
  Broad's GenePattern server.}

  \item{show.row.names}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{rows.to.highlight}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{row.highlight.color}{see HeatMapImage help at
  Broad's GenePattern server.}

  \item{color.scheme}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{color.palette}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{use.color.gradient}{see HeatMapImage help at
  Broad's GenePattern server.}

  \item{memory}{see HeatMapImage help at Broad's
  GenePattern server.}

  \item{verbose}{logical: verbose messages?}

  \item{overwrite}{overwrite existing images?}
}
\value{
  none. creates heatmaps.
}
\description{
  You can either connect to GenePattern server, and submit
  HeatMapImage jobs to the server, or directly execute the
  jar file locally,
}
\details{
  Shout out to Joshua Gould, Broad Institute who wrote the
  jar in the first place.
}
\examples{
# not run
## 1) running locally:
# genepattern.HeatMapImage(in.file=<myfile.gct>, method="local", format="png")
## 2)
# library(GenePattern)
# gp.con <- gp.login("http://pwbc.garvan.unsw.edu.au/gp", "<username>", "<password>")
# genepattern.HeatMapImage(in.file=<myfile.gct>, method="server", gp.connection=gp.con, format="png")
}
\author{
  Mark Cowley, 2009-08-19
}

