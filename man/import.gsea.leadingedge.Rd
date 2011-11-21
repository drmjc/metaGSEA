\name{import.gsea.leadingedge}
\alias{import.gsea.leadingedge}
\title{Import the leading edge subset of genes from each category.
 Must provide the path to the top-level directory that contains the \dQuote{index.html}}
\usage{
  import.gsea.leadingedge(dir, plot = FALSE)
}
\arguments{
  \item{dir}{the top level dir that contains the index.html
  file}

  \item{plot}{logical: leading edge HCL's on top 50 and
  FDR<0.05 genesets? see
  \code{\link{plot.gsea.leadingedge.HCL}}}
}
\value{
  a named list where each element is a vector of gene
  symbols that were contained within the leading edge
  subset for that category. The list names are the gene set
  names
}
\description{
  Import the leading edge subset of genes from each
  category.  Must provide the path to the top-level
  directory that contains the \dQuote{index.html}
}
\author{
  Mark Cowley, 2009-04-06
}

