# How many genesets in a GSEA object?
#
# Parameters:
#	x: a GSEA object
#
# Value:
#	an integer indicating how many genesets are present
#
# Mark Cowley, 2009-10-16


##' How many genesets in a GSEA object?
##' 
##' @param x a GSEA object
##' @return an integer indicating how many genesets are present
##' @author Mark Cowley, 2009-10-16
##' @export
gsea.length <- function(x) {
	nrow(x$tt)
}
