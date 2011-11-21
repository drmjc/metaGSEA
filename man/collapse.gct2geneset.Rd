\name{collapse.gct2geneset}
\alias{collapse.gct2geneset}
\title{Collapse a GCT file from 1 row per gene to 1 row per geneset.}
\usage{
  collapse.gct2geneset(gct.file, gmt.file,
    method = c("mean", "median", "min", "max", "sum"),
    min.size = 1, max.size = 1e+05, gct.out = NULL)
}
\arguments{
  \item{gct.file}{the path to a gct file}

  \item{gmt.file}{the path to a gmt file. untested with a
  gmx file}

  \item{method}{one of mean, median, min, max, sum}

  \item{min.size}{minimum size threshold for each geneset.
  default=1, ie include all genesets with >= 1 gene. GSEA
  default = 15}

  \item{max.size}{maximum size threshold for each geneset.
  default=1e05, ie include all genesets. GSEA default =
  500}

  \item{gct.out}{the path to the output gct file. This will
  overwrite any previous file there.}
}
\value{
  side effect of creating a file at gct.out, with 1 row per
  geneset, as long as the geneset passed the size
  thresholds.  Todo: add PGSEA, maxmean, maxvar as an
  option See also: collapse.gct
}
\description{
  This function collapses gene data to geneset data. NB:
  the gct[,1] must be the same as the ID's in the gmt file.
  Usually, GCT files have probe names in column 1, and GMT
  files have gene symbols. Thus you need to collapse the
  GCT file to 1 row per gene symbol (see collapse.gct).
  Currently supported methods include: mean, median which
  are good if you expect all genes in the pathway to be
  similarly changed min, max if you're looking at pathways
  with high or low abundance sum which is good if the gct
  file is an adjacency matrix (eg summing up the number of
  genes that are mutated in a pathway)
}
\author{
  Mark Cowley, 2011-02-27
}

