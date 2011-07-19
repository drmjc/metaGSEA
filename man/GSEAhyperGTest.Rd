\name{GSEAhyperGTest}
\alias{GSEAhyperGTest}
\title{Perform a Hypergeometric test, upon a GMT file.}
\usage{GSEAhyperGTest(genes, gmt, adjust.method="none")
}
\description{Perform a Hypergeometric test, upon a GMT file.
This is a convenience wrapper around \code{\link{GMThyperGTest}}.}
\value{a \code{data.frame} containing these columns: 
\dQuote{GeneSet}, \dQuote{x}, \dQuote{m}, \dQuote{n}, \dQuote{k}, \dQuote{Pval}, \dQuote{FDR},
where:\cr
\dQuote{x} = the number of white balls drawn without replacement from an urn which contains both black and white balls.\cr
\dQuote{m} = the number of white balls in the urn.\cr
\dQuote{n} = the number of black balls in the urn.\cr
\dQuote{k} = the number of balls drawn from the urn.\cr}
\seealso{\code{\link[stats]{phyper}}}
\author{Mark Cowley, 2009-01-15}
\arguments{\item{genes}{a vector of gene symbols}
\item{gmt}{a gmt object. see \code{\link{import.gsea.gmt}}}
\item{adjust.method}{the multiple hypothesis correction value. see \code{\link{p.adjust}}}
}

