\name{gsea.gct.select.features.rows}
\alias{gsea.gct.select.features.rows}
\title{Subset the rows of a gct file.}
\usage{gsea.gct.select.features.rows(gct.file, gct.out, rows, rows.file)
}
\description{Subset the rows of a gct file.
The GenePattern SelectFeaturesRows tool only operates on the probe id column,
however GCT files also have a gene symbol column.
Restrict the rows of the GCT file, using either a vector or file of probe ID's or gene symbols,
using exact word matches}
\value{none. creates a gct file with rows that match the symbols or probe ID's in 
\code{rows} or \code{rows.file}.}
\author{Mark Cowley, 2010-07-06}
\arguments{\item{gct.file}{the path to a gct file}
\item{gct.out}{the path to the gct result file}
\item{rows}{a character vector of probe id or symbols, or}
\item{rows.file}{a file containing probe id's or symbols. this overrides \code{rows}}
}

