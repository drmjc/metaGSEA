#' Automated GSEA HCL
#' 
#' Automated method for doing a hierarchical clustering of genesets, based on
#' each sets leading edge genes, for all gsea results in a list of gsea
#' results, over a range of different best N, FDR, P, FWER thresholds, and for
#' up/down genesets.
#' 
#' @param x either a gsea list, or a list of gsea lists. see import.gsea
#' @param file either the path to a pdf file which will be created, or NULL to
#'   use the current plotting device
#' @param N vectors of 0 or more thresholds to use. eg N=c(50, 100) means plot
#'   the top 50, or the
#' @param FDR vectors of 0 or more thresholds to use. eg N=c(50, 100) means
#'   plot the top 50, or the
#' @param P vectors of 0 or more thresholds to use. eg N=c(50, 100) means plot
#'   the top 50, or the
#' @param FWER vectors of 0 or more thresholds to use. eg N=c(50, 100) means
#'   plot the top 50, or the top 100 genesets
#' @param max.distance (experimental) threshold to filter out genesets that
#'   are poorly connected to each other. This is a distance threshold, so 1.0
#'   means all genesets, 0.9 is a good place to start.
#' @return makes a number of plots to a pdf file, or the current plotting
#'   device 
#' 
#' @seealso \code{\link{plot_gsea.leadingedge.HCL}}, \code{\link{gsea.leadingedge.distance}}
#' 
#' @author Mark Cowley, 2009-09-03
#' @export
plot_gsea.leadingedge.HCL.auto <- function(x, file=NULL, 
	N=c(50), FDR=c(0.25), P=NULL, FWER=NULL, 
	max.distance=1.0) {

	if( "leading.edge" %in% names(x) )
		x <- list(gsea.collection=x)

	if( !is.null(file) ) {
		pdf.A4(file)
		on.exit(dev.off())
	}

	collections <- names(x)
	for( collection in collections ) {
		cat(".")
		
		main <- "Leading Edge Similarities"
		if( collection != "gsea.collection")
			main <- paste(main, collection, sep=" - ")
			
		# N
		for(n in N) { # skips the loop if N is NULL
			for(direction in c("up", "down")) {
				plot_gsea.leadingedge.HCL(x[[collection]], N=n, direction=direction, main=main, max.distance=max.distance)
			}
		}
		# FDR
		for(fdr in FDR) {
			for(direction in c("up", "down")) {
				plot_gsea.leadingedge.HCL(x[[collection]], FDR=fdr, direction=direction, main=main, max.distance=max.distance)
			}
		}
		# P
		for(p in P) {
			for(direction in c("up", "down")) {
				plot_gsea.leadingedge.HCL(x[[collection]], P=p, direction=direction, main=main, max.distance=max.distance)
			}
		}
		# FWER
		for(fwer in FWER) {
			for(direction in c("up", "down")) {
				plot_gsea.leadingedge.HCL(x[[collection]], FWER=fwer, direction=direction, main=main, max.distance=max.distance)
			}
		}
	}
	cat("\n")
}
