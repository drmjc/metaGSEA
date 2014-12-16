#' Import a gmx file that is commonly used by Broad's GSEA.
#' 
#' Either specify the full path to a gmx file, or from a GSEA PreRanked
#' run, provide the path to the top-level directory that contains the
#' index.html
#' 
#' @param x either the path to a gmx file, the dir that contains index.html
#'   from a GSEA PreRanked run, or a URL to the Broad Institute's website (see
#'   import.gsea.rpt)
#' 
#' @return a named list where each element is a vector of gene symbols the
#'   names are the geneset names
#' 
#' @author Mark Cowley, 2012-12-05
#' @export
import.gsea.gmx <- function(x) {
	if( is.dir(x) ) {
		if( file.exists(file.path(x, "index.html")) ) {
			f <- file.path(x, "edb", "gene_sets.gmx")
			stopifnot( file.exists(f) )
			return( import.gsea.gmx(f) )
		}
		else( stop("Must specify either the gmx file itself, or the dir containing index.html.\n") )
	}
	else if( is.url(x) ) {
		f <- tempfile()
		cat("Attempting to download .gmx/.gmx file.\n")
		download.file(x, f)
		import.gsea.gmx(f)
	}
	else if( is.file(x) ) {
		raw <- read.delim(x, check.names=FALSE, stringsAsFactors=FALSE)
		raw <- raw[-1, ]
		res <- as.list(raw)
		res <- lapply(res, function(x) x[!is.na(x) & nchar(x) > 0])
		names(res) <- colnames(raw)

		return( res )
	}
}
