\name{plot.gsea.venn}
\alias{plot.gsea.venn}
\title{Compare the significant genesets from 2 or 3 GSEA runs using Venn diagrams...}
\usage{plot.gsea.venn(..., gsea.list, names, main="", N, NES, FDR, P, FWER, direction=c("either", "up", "down"))
}
\description{Compare the significant genesets from 2 or 3 GSEA runs using Venn diagrams}
\value{a 2D or 3D venn diagram}
\author{Mark Cowley, 2010-10-15}
\arguments{\item{\dots}{either 2 or 3 GSEA objects by name, or as elements of gsea.list}
\item{gsea.list}{either 2 or 3 GSEA objects by name, or as elements of
gsea.list}
\item{names}{a character vector of names, 1 per GSEA run}
\item{main}{plot title prefix, eg "Venn"}
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

