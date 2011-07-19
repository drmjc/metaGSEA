#' Import the leading edge subset of genes from each category.
#'  Must provide the path to the top-level directory that contains the \dQuote{index.html}
#' 
#' @param dir the top level dir that contains the index.html file
#' @param plot logical: leading edge HCL's on top 50 and FDR<0.05 genesets? see \code{\link{plot.gsea.leadingedge.HCL}}
#' @return a named list where each element is a vector of gene symbols that
#'   were contained within the leading edge subset for that category. The list
#'   names are the gene set names
#' @author Mark Cowley, 2009-04-06
#' @export
import.gsea.leadingedge <- function(dir, plot=FALSE) {
	if( !file.exists(file.path(dir, "index.html")) )
		stop("You did not specify the top level dir that contains index.html")
	
	# import the gmt from this analysis where all gmt elements of rnk
	gmt <- import.gsea.gmt(dir)
	# import the GSEA top table
	tt <- import.gsea.topTable(dir)
	tt <- tt[order(tt$NES, decreasing=TRUE), ]
	# as we iterate through each row of tt, gmt should be in the same order
	gmt <- gmt[match(tt$NAME, names(gmt))]

	# import the rnk, which will let us determine the order of the genes
	rnk <- import.gsea.rnk(dir)
	# reorder the genes within the gmt to match the rnk.

	gmt <- gsea.reorder.gmt.by.rnk(gmt, rnk)

	N <- nrow(tt)
	for(i in 1:N) {
		geneset <- names(gmt)[i]
		# genes <- names(rnk)[ names(rnk) %in% gmt[[geneset]] ]
		genes <- gmt[[geneset]]
		if( tt$DIRECTION[i] =="down" ) {
			genes <- rev(genes)
		}
		if( tt$LEADING.EDGE.SIZE[i] > length(genes) )
			warning(sprintf("%s has fewer genes than the length of the leading edge", geneset))
		gmt[[i]] <- genes[ 1:tt$LEADING.EDGE.SIZE[i] ]
	}
	
	# plot??
	if( plot ) {
		MAIN <- gsea.which.gmt(dir)
		f <- file.path(dir, "LeadingEdge.HCL.pdf")
		pdf.A4(f)
		plot.gsea.leadingedge.HCL(gmt, tt, N=100, main=MAIN)
		plot.gsea.leadingedge.HCL(gmt, tt, FDR=0.05, main=MAIN)
		# plot.gsea.leadingedge.HCL(gmt, tt$pos, FDR=0.1, main=MAIN)
		# plot.gsea.leadingedge.HCL(gmt, tt$pos, FDR=0.25, main=MAIN)
		dev.off()
	}
	
	return( gmt )
}
