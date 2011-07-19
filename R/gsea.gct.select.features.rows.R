#' Subset the rows of a gct file.
#' The GenePattern SelectFeaturesRows tool only operates on the probe id column,
#' however GCT files also have a gene symbol column.
#' Restrict the rows of the GCT file, using either a vector or file of probe ID's or gene symbols,
#' using exact word matches
#' 
#' @param gct.file the path to a gct file
#' @param gct.out the path to the gct result file
#' @param rows a character vector of probe id or symbols, or
#' @param rows.file a file containing probe id's or symbols. this overrides \code{rows}
#' @return none. creates a gct file with rows that match the symbols or probe ID's in 
#' \code{rows} or \code{rows.file}.
#' @author Mark Cowley, 2010-07-06
#' @export
gsea.gct.select.features.rows <- function(gct.file, gct.out, rows=NULL, rows.file=NULL) {
	if( is.null(rows) && is.null(rows.file) )
		stop("Must supply either the rows vector, or rows.file.\n")
	if( is.null(rows) ) {
		rows <- readLines(rows.file)
	}
	gct <- import.gsea.gct(gct.file)
	
	#
	# which column of the gct file should we be searching for the matches?
	#
	overlaps <- list(
		Name=intersect(rows, gct$Name), 
		Description=intersect(rows, gct$Description)
	)
	# which column had the most number of matches to the query?
	column <- which.max(sapply(overlaps,length))[1]
	idx <- which(gct[,column] %in% overlaps[[column]])
	gct <- gct[idx, ]
	
	export.gsea.gct(gct, file=gct.out)
}
