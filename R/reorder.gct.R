#' Function to reorder the rows of a GCT object
#' 
#' take the data in a gct file & reorder the rows by different criteria:\cr
#' - var: sort by variance, from high to low\cr
#' - mean: sort by average abundance, from high to low\cr
#' - median: sort by median abundance, from high to low\cr
#' - max: sort by max abundance, from high to low\cr
#' - sum: sort by sum, from high to low (useful if the GCT is an adjacency
#' matrix/matrix of counts, or scores)\cr
#' 
#' @param gct a gct object
#' @param method see details. Default = 'var'
#' @param reverse if \code{TRUE}, then reverse the default sort order (see details for
#'   default sort order)
#' @author Mark Cowley, 2011-03-16
#' @export
#' @seealso \code{\link{reorder.gct.file}}
reorder.gct <- function(gct, method=c("var", "mean", "median", "max", "sum"), reverse=FALSE) {
	N <- ncol(gct)-2
	
	method <- method[1]
	gct <- switch(method,
		var={
			x <- apply(gct[,3:(N+2)], 1, var)
			gct[order(x, decreasing=TRUE), ]
		},
		mean={
			x <- apply(gct[,3:(N+2)], 1, mean)
			gct[order(x, decreasing=TRUE), ]
		},
		median={
			x <- apply(gct[,3:(N+2)], 1, median)
			gct[order(x, decreasing=TRUE), ]
		},
		sum={
			x <- apply(gct[,3:(N+2)], 1, sum)
			gct[order(x, decreasing=TRUE), ]
		},
		max={
			x <- apply(gct[,3:(N+2)], 1, max)
			gct[order(x, decreasing=TRUE), ]
		},
		stop("Unsupported method")
	)
	
	return(gct)
}

#' Reorder the rows of a GCT file.
#' 
#' see \code{\link{reorder.gct}} for more details
#' 
#' @param gct.file the path to a gct file
#' @param gct.file.out the name of the GCT file to be created with the
#'   reordered rows
#' @inheritParams reorder.gct
#' @return creates a file at gct.file.out, where the rows have been reordered
#' @author Mark Cowley, 2011-03-16
#' @export
#' @seealso \code{\link{reorder.gct}}
reorder.gct.file <- function(gct.file, gct.file.out, method=c("var", "mean", "median", "max", "sum"), reverse=FALSE) {
	gct <- import.gsea.gct(gct.file)
	gct <- reorder.gct(gct, method=method, reverse=reverse)
	export.gsea.gct(gct, file=gct.file.out)
}
