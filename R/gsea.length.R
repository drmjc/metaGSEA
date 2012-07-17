#' How many genesets in a GSEA object?
#' 
#' @param x a GSEA object
#' @return an integer indicating how many genesets are present
#' @author Mark Cowley, 2009-10-16
#' @export
gsea.length <- function(x) {
	nrow(x$tt)
}
