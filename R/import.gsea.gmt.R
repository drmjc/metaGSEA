#' Import a gmt file that is commonly used by Broad's GSEA.
#' 
#' Either specify the full path to a gmt file, or from a GSEA PreRanked
#' run, provide the path to the top-level directory that contains the
#' index.html
#' 
#' @param x either the path to a gmt file, the dir that contains index.html
#'   from a GSEA PreRanked run, or a URL to the Broad Institute's website (see
#'   import.gsea.rpt)
#' 
#' @return a named list where each element is a vector of gene symbols the
#'   names are the geneset names
#' 
#' @author Mark Cowley, 2009-01-14
#' @export
import.gsea.gmt <- function(x) {
	if( is.dir(x) ) {
		if( file.exists(file.path(x, "index.html")) ) {
			f <- file.path(x, "edb", "gene_sets.gmt")
			stopifnot( file.exists(f) )
			return( import.gsea.gmt(f) )
		}
		else( stop("Must specify either the gmt file itself, or the dir containing index.html.\n") )
	}
	else if( is.url(x) ) {
		f <- tempfile()
		cat("Attempting to download .gmt/.gmx file.\n")
		download.file(x, f)
		import.gsea.gmt(f)
	}
	else if( is.file(x) ) {
		raw <- readLines(x, warn=FALSE) # !warn so that missing final EOL is ignored
		res <- strsplit(raw, "\t")
		names <- sapply(res, "[", 1)
		res <- lapply(res, function(x) x[3:length(x)])
		names(res) <- names

		return( res )
	}
}
