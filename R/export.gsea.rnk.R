#' Export a vector of values as a GSEA pre-ranked list.
#' 
#' @param values a numeric vector
#' @param names the names associated with each value. If \code{NULL}, then the
#'   \code{names(values)} are used
#' @param file the output file name. No hyphens, or a warning will be issued.
#' @param resort logical: re-sort the values such that they are decreasing? see order(,
#'   decreasing=TRUE)
#' @param sigfigs the number of decimal places, defaults to 4.
#' @return <none>, makes a file
#' @author Mark Cowley, 2008-09-08
#' @export
export.gsea.rnk <- function(values, names=NULL, file=NULL, resort=TRUE, sigfigs=4) {
	if( grepl("-", file) ) warning("file contains hyphens, which GSEA will probably complain about.")

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
