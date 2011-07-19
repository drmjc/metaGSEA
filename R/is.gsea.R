# Is x a GSEA object
#
# Parameters:
#	x: an object. NOT a list of of objects
#
# Value:
#	logical.
#
# Mark Cowley, 2009-12-16


##' Is x a GSEA object
##' 
##' @param x an object. NOT a list of of objects
##' @return logical.
##' @author Mark Cowley, 2009-12-16
##' @export
is.gsea <- function(x) {
	all( c("tt", "leading.edge", "rpt") %in% names(x) )
}
