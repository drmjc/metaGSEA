#' Export a GSEA rpt file
#' 
#' @param rpt a rpt object (ie a list). see import.gsea.rpt
#' @param f the filename
#' @return none. creates a rpt file
#' @author Mark Cowley, 2009-12-11
#' @export
#' 
export.gsea.rpt <- function(rpt, f) {
	if( ! "producer_class" %in% names(rpt) )
		stop("argument does not look like a gsea rpt object.\n")

	res <- t(list2df(rpt))
	res <- rownames2col(res, 1, "name")
	res <- cbind("param", res)
	res <- rowPaste(res, sep="\t")
	res <- sub("param\tproducer", "producer", res)
	res <- sub("param\tfile", "file", res)
	o <- c(grep("^producer", res), grep("^param", res), grep("^file", res))
	res <- res[o]
	
	writeLines(res, f)
}
