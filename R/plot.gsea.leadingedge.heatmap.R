#' Generate a heatmap of genesets from a GSEA run
#' 
#' Inspired by the heatmap plot in the BROAD's GSEA Leading Edge Viewer tool,
#' this function creates a plot
#' which represents the similarity/dis-similarity between genesets with
#' colours. The colour is based on the
#' \code{Jacquard}. see \code{\link{dist.gsea}}
#' @note
#' Recommend saving to file using pdf(..., width=12, height=12)
#' 
#' @param x a GSEA object
#' @param main the plot title
#' @param gridcolour the colour of the grid to overlay on top of the image. 
#'     see sepcolor parameter of \code{\link[gplots]{heatmap.2}}
#' @return nothing. generates a plot.
#' @author Mark Cowley, 2009-10-29
#' @export
plot.gsea.leadingedge.heatmap <- function(x, main="Geneset similarity", gridcolour="lightgrey") {
	tmp <- dist.gsea(x$leading.edge)
	tmp <- as.matrix(tmp)
	tmp <- 1-tmp
	tmp[upper.tri(tmp, diag=FALSE)] <- 0
	# upperTriangle(tmp) <- 0 # remove dependency on gdata for such a simple function

	# image.table(tmp, col=colour.step("white", "green", steps=21), zlim=c(0,1), grid.col="white", ylabels=colnames(tmp), xlabels=colnames(tmp), legend=FALSE, xlab="", ylab="", main="", xlabels.pos="top", ylabels.pos="right")

	require(gplots)

	# pdf(file.pdf, 12, 12)
	heatmap.2(tmp, Rowv=FALSE, Colv=FALSE, dendrogram="none", scale="none", col=colour.step("white", "green", steps=21), trace="none", margins=c(25, 25), density.info="histogram", keysize=0.9, 
	rowsep=1:ncol(tmp), colsep=1:nrow(tmp), sepwidth=c(0.005,0.005), sepcolor=gridcolour,
	main=main)
	# dev.off()
}

#' Generate a heatmap of genesets from a GSEA run, across a number of standard thresholds
#' The default settings will generate 4 heatmaps, of the up/down, top 50, and FDR<0.25 genesets.
#' 
#' @param x a GSEA obect
#' @param file the path to a PDF file, or \code{NULL} to use \code{dev.cur}
#' @param N an integer vector of thresholds on the top N genesets, sorted by |NES|
#' @param FDR an integer vector of FDR thresholds
#' @param P an integer vector of unadjusted P-value thresholds
#' @param FWER an integer vector of FWER thresholds
#' @param pdf.width the width of the pdf device in inches
#' @param pdf.height the height of the pdf device in inches
#' @return nothing. creates a pdf file.
#' @author Mark Cowley, 2010-10-14
#' @export
plot.gsea.leadingedge.heatmap.auto <- function(x, file=NULL, 
	N=c(50), FDR=c(0.25), P=NULL, FWER=NULL, pdf.width=12, pdf.height=12) {

	if( "leading.edge" %in% names(x) )
		x <- list(gsea.collection=x)

	if( !is.null(file) ) {
		pdf(file, width=pdf.width, height=pdf.height)
		on.exit(dev.off())
	}

	collections <- names(x)
	for( collection in collections ) {
		cat(".")
		
		main <- "Geneset Similarities"
		if( collection != "gsea.collection")
			main <- paste(main, collection, sep=" - ")
			
		# N
		for(n in N) { # skips the loop if N is NULL
			for(direction in c("up", "down")) {
				tmp.main <- paste(main, " - ", direction, " - topN=", n, sep="")
				tmp <- gsea.filter(x[[collection]], N=n, direction=direction)
				if( length(tmp$leading.edge) >0 )
					plot.gsea.leadingedge.heatmap(tmp, main=tmp.main)
				else
					plot.blank(main=tmp.main, box=TRUE, message="**** too few genesets pass thresholds ****")
			}
		}
		# FDR
		for(fdr in FDR) {
			for(direction in c("up", "down")) {
				tmp.main <- paste(main, " - ", direction, " - FDR<", fdr, sep="")
				tmp <- gsea.filter(x[[collection]], FDR=fdr, direction=direction)
				if( length(tmp$leading.edge) >0 )
					plot.gsea.leadingedge.heatmap(tmp, main=tmp.main)
				else
					plot.blank(main=tmp.main, box=TRUE, message="**** too few genesets pass thresholds ****")
			}
		}
		# P
		for(p in P) {
			for(direction in c("up", "down")) {
				tmp.main <- paste(main, " - ", direction, " - P<", p, sep="")
				tmp <- gsea.filter(x[[collection]], P=p, direction=direction)
				if( length(tmp$leading.edge) >0 )
					plot.gsea.leadingedge.heatmap(tmp, main=tmp.main)
				else
					plot.blank(main=tmp.main, box=TRUE, message="**** too few genesets pass thresholds ****")
			}
		}
		# FWER
		for(fwer in FWER) {
			for(direction in c("up", "down")) {
				tmp.main <- paste(main, " - ", direction, " - FWER<", fwer, sep="")
				tmp <- gsea.filter(x[[collection]], FWER=fwer, direction=direction)
				if( length(tmp$leading.edge) >0 )
					plot.gsea.leadingedge.heatmap(tmp, main=tmp.main)
				else
					plot.blank(main=tmp.main, box=TRUE, message="**** too few genesets pass thresholds ****")
			}
		}
	}
	cat("\n")
}
