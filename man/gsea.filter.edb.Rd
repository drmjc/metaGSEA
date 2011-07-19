\name{gsea.filter.edb}
\alias{gsea.filter.edb}
\title{Filter a GSEA edb file, restricting the entries to restricted set of...}
\usage{gsea.filter.edb(x, genesets)
}
\description{Filter a GSEA edb file, restricting the entries to restricted set of
genesets.}
\value{a gsea list with the edb object restricted to only those that were
specified by \code{genesets}}
\author{Mark Cowley, 2009-10-13}
\arguments{\item{x}{a gsea list (not just the edb file)}
\item{genesets}{a vector of geneset names. Minimal error checking done, so
make sure they exist in your edb file!}
}

