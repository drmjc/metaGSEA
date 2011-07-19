\name{gsea.rename.genesets}
\alias{gsea.rename.genesets}
\title{Rename the genesets from within a GSEA object.}
\usage{gsea.rename.genesets(x, prefix, suffix, fdr=FALSE, rank=FALSE, direction=FALSE, maxlen=60)
}
\description{Rename the genesets from within a GSEA object.}
\details{Rename genesets using a prefix and suffix, such that COWLEY_UP becomes:
PREFIX_COWLEY_UP_SUFFIX
Optionally, the rank, FDR and direction can also be included so that
COWLEY_UP becomes:
1: PREFIX_COWLEY_UP_SUFFIX UP (0.002)}
\value{a GSEA list with modified names.}
\author{Mark Cowley, 2009-10-14}
\arguments{\item{x}{a GSEA list, optionally with an edb. see import.gsea(...,
edb=TRUE)}
\item{prefix}{an optional prefix (such as the treatment, eg WTvsTreatA).}
\item{suffix}{an optional prefix (such as the treatment, eg WTvsTreatA).}
\item{fdr}{add the fdr in parentheses to the end? eg " (0.001)"}
\item{rank}{add the rank at the start? eg "12: "}
\item{direction}{add UP or DOWN to the end?}
\item{maxlen}{the max number of characters in the geneset name (before we
add rank & fdr & direction to it)}
}

