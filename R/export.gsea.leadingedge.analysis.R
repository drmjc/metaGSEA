#' Export the plots and tables required for a GSEA leading edge analysis.
#' 
#' @param x a GSEA object. Usually this will have been heavily subset to
#'   genesets of interest. eg using gsea.split, or gsea.filter(...,
#'   gensets=<>)
#' @param dir the directory where files will be made.
#' @param prefix an optional prefix to differentiate these files from others
#'   that you might make. eg WTvsKO1
#' @param sep the seperator to use in constructing the outputfilename. 
#'    Ignored if \code{prefix=NULL}
#' @return Creates N files:\cr
#' leadingedge.plots.pdf: 4 plots in a single PDF file\cr
#' leadingedge.genes.txt: simple text format of genesymbols. useful for
#'   GO analysis (using eg DAVID)\cr
#' leadingedge.genecounts.xls: a 2 column table
#'   of genes with their counts in the N genesets\cr
#' leadingedge.adjmat.xls: a
#'   genes x genesets adjacency matrix.\cr
#' leadingedge.tree.cdt: a fileformat
#'   that lets you view the hierarchical clustering of genesets using the
#'   HierarchicalClusteringViewer tool @@ Broad's GenePattern.\cr
#'   leadingedge.tree.gtr: used in conjunction with the cdt file.
#' @author Mark Cowley, 2009-10-29
#' @export
export.gsea.leadingedge.analysis <- function(x, dir, prefix=NULL, sep="-") {
	if( is.null(prefix) ) {
		prefix <- ""
		sep <- ""
	}
	
	# export the plots
	f <- file.path(dir, paste(prefix, "leadingedge.plots.pdf", sep=sep))
	pdf(f, 12, 12)
	plot_gsea.leadingedge(x, main=prefix)
	dev.off()
	
	# export a list of non-unique genesymbols (for GO analysis)
	genes <- sort(unique(unlist(x$leading.edge)))
	f <- file.path(dir, paste(prefix, "leadingedge.genes.txt", sep=sep))
	writeLines(genes, f)
	
	# export gene counts
	uc <- ucounts(unlist(x$leading.edge))
	uc <- uc[order(uc, decreasing=TRUE)]
	uc <- cbind(GeneSymbol=names(uc), Count=uc)
	f <- file.path(dir, paste(prefix, "leadingedge.genecounts.xls", sep=sep))
	write.xls(uc, f, row.names=FALSE, col.names=TRUE)
	
	# export adjmat
	am <- gsea2adjmat(x)
	f <- file.path(dir, paste(prefix, "leadingedge.adjmat.xls", sep=sep))
	write.xls(am, f,  row.names="GeneSymbol", col.names=TRUE)
	
	# export the cdt/gtr for viewing the HCL in genepattern.
	f <- file.path(dir, paste(prefix, "leadingedge.tree.cdt", sep=sep))
	export.gsea.cdt(x, f, gtr=TRUE, rename=FALSE)	
	
	# export the combined top table as a real XLS
	f <- file.path(dir, paste(prefix, "gsea_report.xls", sep=sep))
	write.xls(x$tt, f)

	# export a GSEA summary
	f <- file.path(dir, paste(prefix, "gsea_summary.xls", sep=sep))
	write.xls(gsea.summarise(x), f, row.names=TRUE, col.names=TRUE)
	
}
