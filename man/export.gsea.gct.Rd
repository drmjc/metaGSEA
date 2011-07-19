\name{export.gsea.gct}
\alias{export.gsea.gct}
\title{Function to export a GenePattern gct file.}
\usage{export.gsea.gct(data, description, file, chip, round=4, version="#1.2", missing="", ...)
}
\description{Function to export a GenePattern gct file.}
\details{\code{data} can be either a \code{matrix} or \code{data.frame} of numbers, 
in which case description and or chip should be provided to annotate the rows; 
or a \code{data.frame} containing a \sQuote{Name} 
and \dQuote{Description} columns (see \code{\link{import.gsea.gct}})}
\author{Mark Cowley, 2008-08-07}
\arguments{\item{data}{a data.frame WITH rownames, and optionally can be from
import.gsea.gct; ie first 2 columns are Name, Description}
\item{description}{a vector of annotations. IF it has names, then these
names must match the rownames of data, in which case we make sure they're
in the same order. If no names, then we assume that they're in the same
order. If NULL, then you must either provide a chip, or have a
Description column in the data.}
\item{file}{the output file name}
\item{chip}{instead of specifying description, you can specify a GSEA chip
object, and a description will be made for you. This overrides
description. Default is NULL to ignore.}
\item{round}{the number of digits to round the numbers to - defaults to 3}
\item{version}{The first line of the gct file will have this tag.}
\item{missing}{the string to use for missing data (in the expression data).}
\item{\dots}{Currently unused.}
}

