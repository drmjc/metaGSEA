#' Determine the distance between various leading edge genes, in terms of the
#' jacquard.
#' 
#' @param x a list of vectors of gene symbols
#' @return an object of class \sQuote{dist} where each value is \code{1-Jacquard(a, b)}
#'   where \code{a} and \code{b} are the leading edges from two genesets.
#' @author Mark Cowley, 2009-04-06
#' @seealso \code{\link{dist}}, \code{\link{dist.gsea}}
#' @export
gsea.leadingedge.distance <- function(x) {
	res <- matrix(0, length(x), length(x))
	for(row in 1:length(x)) {
		for(column in 1:row) {
			res[row,column] <- 1 - jacquard(x[[row]], x[[column]])
		}
	}
	if( is.null(names(x)) ) {
		rownames(res) <- 1:length(x)
		colnames(res) <- 1:length(x)
	}
	else {
		rownames(res) <- colnames(res) <- names(x)
	}
	# diag(res) <- 0
	# res <- res[-1, -ncol(res)]
	res <- as.dist(res)
	
	return( res )
}

#' Determine the distance between various leading edge genes, in terms of the
#' jacquard.
#' 
#' This is a convenience wrapper around \code{\link{gsea.leadingedge.distance}}
#' 
#' @param x a list of vectors of gene symbols
#' @return an object of class 'dist' where each value is 1-Jacquard(a, b)
#'   where a and b are the leading edges from two genesets.
#' @author Mark Cowley, 2009-04-06
#' @seealso \code{\link{dist}}, \code{\link{gsea.leadingedge.distance}}
#' @export
dist.gsea <- function(x) {
	gsea.leadingedge.distance(x)
}
