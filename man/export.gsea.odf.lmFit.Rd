\name{export.gsea.odf.lmFit}
\alias{export.gsea.odf.lmFit}
\title{Export a limma lmFit to a GSEA ODF file.}
\usage{export.gsea.odf.lmFit(fit1, fit2, coef=1, file, gct.file, cls.file, description, collapse=FALSE, probe2gene)
}
\description{Export a limma lmFit to a GSEA ODF file.
From two linear model files (see details), the coefficient of interest, export an
ODF file, which can be viewed in the GenePattern ComparativeMarkerSelectionViewer 
tool. The ODF file that is created is heavily modelled on the ODF file that the
ComparativeMarkerSelection GenePattern module creates.}
\details{This is a key function for the LimmaGP GenePattern module. 
In LimmaGP, we always perform 
two linear model fits, the first is a treatment means 
parameterisation (\code{fit1}, hereafter), 
the second is a contrast fit after constructing the appropriate contrasts (\code{fit2}, hereafter).
From \code{fit1}, we can get the group means, and standard errors, and from \code{fit2}, 
we can get the logFC,
moderated t-statistics, pvalues, and FDR/FWER values.}
\value{nothing. writes an ODF file.}
\author{Mark Cowley, 2011-07-19}
\arguments{\item{fit1}{a linear model fit, \emph{before} fitting contrasts, and running eBayes}
\item{fit2}{a linear model fit, \emph{after} fitting contrasts, and running eBayes}
\item{coef}{the coefficient to be exported. see \code{\link[limma]{topTable}}}
\item{file}{the path to the ODF file}
\item{gct.file}{the name of the gct file that the linear models were fitted to. 
This is written to the ODF header}
\item{cls.file}{the name of the cls file that the linear models were fitted to. 
This is written to the ODF header}
\item{description}{a character vector, which will end up in the description column.
if \code{NULL}, then the probe ID's will be used.}
\item{collapse}{logical: whether to collapse probes to genes? if \code{TRUE}, then 
\code{probe2gene} must be set.}
\item{probe2gene}{if \code{collapse=TRUE}, then a 2 column data.frame mapping from
probe to gene symbol should be supplied. default: \code{NULL}.}
}

