# Function to reorder the rows of a GCT file
#
# Details:
#	take the data in a gct file & reorder the rows by different criteria:
#	- var: sort by variance, from high to low
#	- mean: sort by average abundance, from high to low
#	- median: sort by median abundance, from high to low
#	- sum: sort by sum, from high to low (useful if the GCT is an adjacency matrix/matrix of counts, or scores)
#
# Parameters:
#	gct: a gct object
#	method: see details. Defualt = 'var'
#	reverse: if TRUE, the reverse the default sort order (see details for default sort order)
#
# Mark Cowley, 2011-03-16


##' Function to reorder the rows of a GCT file
##' 
##' take the data in a gct file & reorder the rows by different criteria:
##' - var: sort by variance, from high to low
##' - mean: sort by average abundance, from high to low
##' - median: sort by median abundance, from high to low
##' - sum: sort by sum, from high to low (useful if the GCT is an adjacency
##' matrix/matrix of counts, or scores)
##' 
##' @param gct a gct object
##' @param method see details. Defualt = 'var'
##' @param reverse if TRUE, the reverse the default sort order (see details for
##'   default sort order)
##' @author Mark Cowley, 2011-03-16
##' @export
reorder.gct <- function(gct, method=c("var", "mean", "median", "sum"), reverse=FALSE) {
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
		stop("Unsupported method")
	)
	
	return(gct)
}

# Reorder the rows of a GCT file.
#
# Details: see reorder.gct
#
# Parameters:
#	gct.file, method, reverse: see reorder.gct
#	gct.file.out: the name of the GCT file to be created with the reordered rows
#
# Value:
#	creates a file at gct.file.out, where the rows have been reordered
#
# Mark Cowley, 2011-03-16


##' Reorder the rows of a GCT file.
##' 
##' Details: see reorder.gct
##' 
##' @param gct.file see reorder.gct
##' @param method see reorder.gct
##' @param reverse see reorder.gct
##' @param gct.file.out the name of the GCT file to be created with the
##'   reordered rows
##' @return creates a file at gct.file.out, where the rows have been reordered
##' @author Mark Cowley, 2011-03-16
##' @export
reorder.gct.file <- function(gct.file, gct.file.out, method=c("var", "mean", "median", "sum"), reverse=FALSE) {
	gct <- import.gsea.gct(gct.file)
	gct <- reorder.gct(gct, method=method, reverse=reverse)
	export.gsea.gct(gct, file=gct.file.out)
}
