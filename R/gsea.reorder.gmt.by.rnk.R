# Reorder the genes within a gmt object to match the order of genes within a rnk object
#
# Parameters:
#	gmt: a gmt object (ie a named list of genes)
#	rnk: a rnk object (ie named list of gene scores, where genes are the names)
#
# Value:
#	a gmt just like input, but with the genes reordered to match the order inside the rnk
#
# Mark Cowley, 2009-12-07
#


##' Reorder the genes within a gmt object to match the order of genes within a
##' rnk object
##' 
##' @param gmt a gmt object (ie a named list of genes)
##' @param rnk a rnk object (ie named list of gene scores, where genes are the
##'   names)
##' @return a gmt just like input, but with the genes reordered to match the
##'   order inside the rnk
##' @author Mark Cowley, 2009-12-07
##' @export
gsea.reorder.gmt.by.rnk <- function(gmt, rnk) {
	N <- length(gmt)
	rnk <- rnk[order(rnk, decreasing=TRUE)]
	rnk <- names(rnk)
	
	for(i in 1:N) {
		gmt[[i]] <- rnk[ rnk %in% gmt[[i]] ]
	}

	return( gmt )
}
