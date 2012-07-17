#' Export a vector of values as a GSEA pre-ranked list.
#' 
#' @details Specifying up-/down-regulated genes. One of the
#' biggest causes of confusion, especially when revisiting
#' results generated a long time ago, are my positive genes
#' up in class A, or class B? If you specify \code{up.class}
#' and \code{down.class}, a header is written into in the rnk file.
#' The new GseaPreRanked GenePattern module will read that
#' header & label the geneset reports with these names, rather
#' than na_pos and na_neg.
#' 
#' @param values a numeric vector
#' @param names the names associated with each value. If \code{NULL}, then the
#'   \code{names(values)} are used
#' @param file the output file name. No hyphens, or a warning will be issued.
#' @param resort logical: re-sort the values such that they are decreasing? see order(,
#'   decreasing=TRUE)
#' @param sigfigs the number of significant figures, defaults to 4.
#' @param up.class The up-regulated genes are up in \code{up.class}
#' @param down.class The down-regulated genes are down in \code{down.class}
#' 
#' @return <none>, makes a file
#' 
#' @author Mark Cowley, 2008-09-08
#' @export
export.gsea.rnk <- function(values, names=NULL, file=NULL, resort=TRUE, sigfigs=4, up.class=NA, down.class=NA) {
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
	
	# values <- round(values, sigfigs)
	values <- signif(values, sigfigs)
	
	tmp <- data.frame(names, values, stringsAsFactors=FALSE)
	
	f <- file(file, "w")
	if( !is.na(up.class) && !is.na(down.class) ) {
		writeLines(sprintf("#up.class=%2", up.class), f)
		writeLines(sprintf("#down.class=%2", down.class), f)
	}
	write.delim(tmp, f, col.names=FALSE)
	close(f)
}
# CHANGELOG
# 2012-07-04
# - added up.class/down.class parameters and header
# - changed from round to signif
# 