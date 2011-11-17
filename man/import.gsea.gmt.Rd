\name{import.gsea.gmt}
\alias{import.gsea.gmt}
\title{Import a gmt file that is commonly used by Broad's GSEA.}
\usage{
  import.gsea.gmt(x)
}
\arguments{
  \item{x}{either the path to a gmt file, the dir that
  contains index.html from a GSEA PreRanked run, or a URL
  to the Broad Institute's website (see import.gsea.rpt)}
}
\value{
  a named list where each element is a vector of gene
  symbols the names are the geneset names
}
\description{
  Either specify the full path to a gmt file, or from a
  GSEA PreRanked run, provide the path to the top-level
  directory that contains the index.html
}
\author{
  Mark Cowley, 2009-01-14
}

