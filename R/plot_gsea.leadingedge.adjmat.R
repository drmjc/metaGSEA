#' Generate a binary adjacency matrix image (like a heatmap with 2 colours).
#' 
#' The occurence of genes in each geneset is represented by blue squares. Like
#' in the Broad's GSEA Leading Edge Viewer tool, genes are in columns, gensets
#' are in rows.
#' 
#' @param x a GSEA object
#' @param cluster if TRUE then do a 2-way hierarchical clustering, else order
#'   genes by decreasing frequency
#' @param main the title
#' @return none. creates a plot.  Todo: introduce colour by looking up stat in
#'   the rnk object. handle the case of lots of genes nicer.
#' @author Mark Cowley, 2009-10-29
#' @export
#' @importFrom gplots heatmap.2
plot_gsea.leadingedge.adjmat <- function(x, main="Leading Edge Adjacency Matrix", cluster=TRUE) {
	adjmat <- gsea2adjmat(x)
	adjmat <- as.matrix(adjmat)
	adjmat <- t(adjmat)
	
	try( {
		heatmap.2(adjmat, 
			Rowv=cluster, Colv=cluster, dendrogram="none", scale="none", 
			col=c("white", "blue"), trace="none", margins=c(5,25), 
			key=FALSE, keysize=0.1, lwd=c(0,10), lhei=c(5,15),
			# rowsep=1:ncol(tmp), colsep=1:nrow(tmp), sepwidth=c(0.005,0.005), sepcolor="white",
			main=main, 
			cexCol = 0.1 + 1/log10(ncol(adjmat)),
			cexRow = 0.1 + 1/log10(nrow(adjmat))
		) },
	silent=TRUE
	)
}
