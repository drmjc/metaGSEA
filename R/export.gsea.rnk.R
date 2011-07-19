#
# Export a vector of values as a GSEA pre-ranked list.
#
# Parameters:
#	values: a numeric vector
#	names: the names associated with each value. If NULL, then the names(values) are used
#	file: the output file name
#	resort: re-sort the values such that they are decreasing. see order(, decreasing=TRUE)
#	sigfigs: the number of decimal places, defaults to 4.
#
# Value:
#	<none>, makes a file
#
# Mark Cowley, 2008-09-08
#


##' Export a vector of values as a GSEA pre-ranked list.
##' 
##' @param values a numeric vector
##' @param names the names associated with each value. If NULL, then the
##'   names(values) are used
##' @param file the output file name
##' @param resort re-sort the values such that they are decreasing. see order(,
##'   decreasing=TRUE)
##' @param sigfigs the number of decimal places, defaults to 4.
##' @return <none>, makes a file
##' @author Mark Cowley, 2008-09-08
##' @export
export.gsea.rnk <- function(values, names=NULL, file=NULL, resort=TRUE, sigfigs=4) {
	dir <- dirname(file)
	if( !file.exists(dir) ) {
		warning("Creating directory:", dir)
		dir.create(dir, recursive=TRUE)
	}
	
	if( is.null(names) ) {
		if( !is.null(names(values)) )
			names <- names(values)
		else
			stop("need to provide names, or name the elements of the values.\n")
	}

	if( resort ) {
		o <- order(values, decreasing=TRUE)
		values <- values[o]
		names <- names[o]
	}
	
	values <- round(values, sigfigs)
	
	tmp <- data.frame(names, values, stringsAsFactors=FALSE)
	
	write.delim(tmp, file, col.names=FALSE)
}
