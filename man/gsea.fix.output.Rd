\name{gsea.fix.output}
\alias{gsea.fix.output}
\title{Fix the output from a GSEA run.}
\usage{gsea.fix.output(dir, genes.gct, fix.index=TRUE, fix.posneg=TRUE, fix.reports=TRUE,
    debug=FALSE, overwrite.heatmap.images=FALSE, verbose=FALSE)
}
\description{Fix the output from a GSEA run.
GSEA and GseaPreRanked both produce rich HTML reports with hyperlinks to 
annotation databaes, images and heatmaps. GseaPreRanked, misses out on some
of the features of the GSEA reports, and both GSEA and GseaPreRanked output
can be improved somewhat (IMHO).}
\details{Here's a list of the fixes, all to the HTML files:\cr
GSEA and GseaPreRanked:\cr
- Within the GeneSet column, in the geneset reports, remove the link to Affymetrix, and add hyperlinks to
Entrez and GeneCards\cr
- Within the geneset reports, add a hyperlink to the genesets in the header section\cr
- Within the up/down-regulated summaries, add hyperlinks to MSigDB for all genesets, not
just the top N that had individual reports made\cr
- Within the geneset reports, add a gene decscription into the appropriate column, by searching
within the chip file listed in the rpt file\cr


GseaPreRanked only:\cr
- Within \sQuote{index.html}, instead of naming the phenotypes na_pos and na_neg, call 
them up- and down-regulated genesets\cr
- within the individual geneset reports, add a heatmap (if a GCT file is provided)\cr

Control over which elements are fixed:\cr
\code{fix.posneg}: fix the 2 html files that summarise the stats for the
up/down regulated genesets. This only makes sense for
GseaPreRanked output, but no harm in leaving it \code{TRUE} for gsea output.\cr
\code{fix.index}: fix the index.html file? This only makes sense for
GseaPreRanked output, but no harm in leaving it \code{TRUE} for gsea output.\cr
\code{fix.reports}: fix the individual geneset reports? This controls whether the
geneset report fixes reported above, including adding heatmaps is done.\cr

You can import individual, or multiple GSEA runs via the \code{dir} parameter:\cr
- a single gsea directory, eg "./GSEA/pos6_versus_neg6/c1_all.Gsea.1252052322484"\cr
- the parent dir containing multiple GSEA runs, eg "./GSEA/pos6_versus_neg6" which will fix all GSEA dirs within it.\cr
- a single geneset report file, eg "./GSEA/pos6_versus_neg6/c1_all.Gsea.1252052322484/CHR1P21.html"\cr

If you ran GseaPreRanked, then you don't get any
heatmaps. Specifying a gct file via the \code{genes.gct} parameter will
make heatmaps for each of the individual geneset reports. 
NB: the gct file should have one row per
gene symbol only, and the row name should be the gene symbol, NOT the probesetID
TODO: allow a probe-level gct file which matches the probe-level rnk file & which 
is generally easier to obtain a priori.}
\value{none. it just fixes gsea outputs in situ}
\author{Mark Cowley, 2009-07-09}
\arguments{\item{dir}{see details}
\item{genes.gct}{The path to a GCT file. see details.}
\item{fix.index}{logical: fix the index.html file? see details.}
\item{fix.posneg}{logical: Fix the pair of up-/down-regulated geneset summary tables? see details.}
\item{fix.reports}{logical: fix the individual geneset reports? see details}
\item{debug}{logical: print debugging messages?}
\item{overwrite.heatmap.images}{logical: overwrite already existing heatmap images?}
\item{verbose}{logical: print verbose messages?}
}

