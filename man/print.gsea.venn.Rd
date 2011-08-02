\name{print.gsea.venn}
\alias{print.gsea.venn}
\title{Print a Venn diagram of overlaps between genesets from >= 2 GSEA runs...}
\usage{print.gsea.venn(..., gsea.list, names, N, NES, FDR, P, FWER, direction=c("either",
    "up", "down"))
}
\description{Print a Venn diagram of overlaps between genesets from >= 2 GSEA runs}
\value{none. prints a venn diagram}
\author{Mark Cowley}
\arguments{\item{\dots}{>=2 GSEA results; or}
\item{gsea.list}{a list of GSEA objects}
\item{names}{the GSEA run names for labelling the Venn diagram}
\item{N}{see gsea.filter. It doesn't make sense to do a Venn diagram on all
genesets, so you need to filter to find the significant ones.}
\item{NES}{see gsea.filter. It doesn't make sense to do a Venn diagram on
all genesets, so you need to filter to find the significant ones.}
\item{FDR}{see gsea.filter. It doesn't make sense to do a Venn diagram on
all genesets, so you need to filter to find the significant ones.}
\item{P}{see gsea.filter. It doesn't make sense to do a Venn diagram on all
genesets, so you need to filter to find the significant ones.}
\item{FWER}{see gsea.filter. It doesn't make sense to do a Venn diagram on
all genesets, so you need to filter to find the significant ones.}
\item{direction}{see gsea.filter. It doesn't make sense to do a Venn
diagram on all genesets, so you need to filter to find the significant
ones.}
}

