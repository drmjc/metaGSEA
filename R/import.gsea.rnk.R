#' Import a rnk file used for running GSEA PreRanked.
#' 
#' Either specify the full path to a rnk file, or from a GSEA PreRanked
#' run, provide the path to the top-level directory that contains the
#' index.html.
#' 
#' @note this imports the collapsed rnk file.
#' 
#' @param x the path to either a rnk file, or the dir that contains index.html
#' 
#' @return a named numeric vector of values
#' 
#' @author Mark Cowley, 2009-04-28
#' @export
import.gsea.rnk <- function(x) {
	if( is.gsea.dir(x) ) {
		# f <- dir(file.path(x, "edb"), pattern=".*rnk$", full.names=TRUE)
		# stopifnot( file.exists(f) )
		rpt <- import.gsea.rpt(x)
		f <- rpt$collapsed_rnk
		if( !file.exists(f) )
			stop("Can't find the collapsed rnk file, within the edb directory. Have you moved this GSEA directory, and not updated the rpt file?")
		return( import.gsea.rnk(f) )
	}
	else if( is.file(x) && file.exists(x) ) {
		rnk <- read.delim(x, header=FALSE)
		rnk <- rnk[order(rnk[,2], decreasing=TRUE), ]
		res <- rnk[,2]
		names(res) <- rnk[,1]

		return( res )
	}
	else
		stop("Must specify either the rnk file itself, or the dir containing index.html.\n")
}
