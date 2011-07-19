#' Export a GSEA chip file
#' Export a GSEA chip file -- this basically ensures that the correct column
#' names are exported.
#' 
#' Export a GSEA chip file -- this basically ensures that the correct column
#' names are exported.
#' 
#' Assumes that the user provided these columns in order: 1. Probe ID 2. Gene
#' Symbol 3. Gene Description 4. optional Gene Symbol Aliases
#' 
#' @param x a matrix-like object. see details, and
#'   \code{\link{import.gsea.chip}}
#' @param f the output file name
#' @param resort.rows re-order the rows using \code{order} on column 1.
#' @param collapse.rows ensure that each probe is only outputted once
#' @param na.strings what character value to replace all of the NA's with. eg
#'   \dQuote{} or \dQuote{NULL} or \dQuote{---}
#' @return none - creates the file.
#' @author Mark Cowley, 2008-09-01
#' @seealso \code{\link{import.gsea.chip}}
#' @references
#'   \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#chip}
#' @keywords file IO
#' @export
export.gsea.chip <- function(x, f, resort.rows=TRUE, collapse.rows=TRUE, na.strings="---") {
	x[x==""] <- NA
	
	if( resort.rows )
		x <- x[order(x[,1]), ] # resort the rows by probe name
	if( collapse.rows ) 
		x <- x[match(unique(x[,1]), x[,1]), ]
	if( any(is.na(x)) )
		x[is.na(x)] <- na.strings

	OUT <- file(f, "w")
	hdr <- c("Probe Set ID", "Gene Symbol", "Gene Title", "Aliases")[1:ncol(x)]
	write(paste(hdr, collapse="\t"), OUT)
	write.delim(x, OUT, col.names=FALSE, row.names=FALSE)

	close(OUT)
}
