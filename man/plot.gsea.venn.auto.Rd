\name{plot.gsea.venn.auto}
\alias{plot.gsea.venn.auto}
\title{Plot a Venn diagram of overlaps between genesets from >= 2 GSEA runs across...}
\usage{plot.gsea.venn.auto(..., gsea.list, names, file, main="", N=c(50), FDR=c(0.25), P, FWER)
}
\description{Plot a Venn diagram of overlaps between genesets from >= 2 GSEA runs across
a range of thresholds.
By default, this will plot 4 venn diagrams: top 50 up, top 50 down, FDR<0.25 up,
FDR<0.25 down}
\value{none. prints a venn diagram}
\author{Mark Cowley}
\arguments{\item{\dots}{>=2 GSEA results; or}
\item{gsea.list}{a list of GSEA objects}
\item{names}{the GSEA run names for labelling the Venn diagram}
\item{file}{the path to the output pdf file, or \code{NULL}.}
\item{main}{the plot title}
\item{N}{see gsea.filter. It doesn't make sense to do a Venn diagram on all
genesets, so you need to filter to find the significant ones.}
\item{FDR}{see gsea.filter. It doesn't make sense to do a Venn diagram on
all genesets, so you need to filter to find the significant ones.}
\item{P}{see gsea.filter. It doesn't make sense to do a Venn diagram on all
genesets, so you need to filter to find the significant ones.}
\item{FWER}{see gsea.filter. It doesn't make sense to do a Venn diagram on
all genesets, so you need to filter to find the significant ones.}
}

