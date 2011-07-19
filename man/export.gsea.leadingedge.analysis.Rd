\name{export.gsea.leadingedge.analysis}
\alias{export.gsea.leadingedge.analysis}
\title{Export the plots and tables required for a GSEA leading edge analysis.}
\usage{export.gsea.leadingedge.analysis(x, dir, prefix, sep="-")
}
\description{Export the plots and tables required for a GSEA leading edge analysis.}
\value{Creates N files:\cr
leadingedge.plots.pdf: 4 plots in a single PDF file\cr
leadingedge.genes.txt: simple text format of genesymbols. useful for
GO analysis (using eg DAVID)\cr
leadingedge.genecounts.xls: a 2 column table
of genes with their counts in the N genesets\cr
leadingedge.adjmat.xls: a
genes x genesets adjacency matrix.\cr
leadingedge.tree.cdt: a fileformat
that lets you view the hierarchical clustering of genesets using the
HierarchicalClusteringViewer tool @ Broad's GenePattern.\cr
leadingedge.tree.gtr: used in conjunction with the cdt file.}
\author{Mark Cowley, 2009-10-29}
\arguments{\item{x}{a GSEA object. Usually this will have been heavily subset to
genesets of interest. eg using gsea.split, or gsea.filter(...,
gensets=<>)}
\item{dir}{the directory where files will be made.}
\item{prefix}{an optional prefix to differentiate these files from others
that you might make. eg WTvsKO1}
\item{sep}{the seperator to use in constructing the outputfilename. 
Ignored if \code{prefix=NULL}}
}

