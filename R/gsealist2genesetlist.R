#' Extract the geneset names from a list of GSEA objects.
#' 
#' From a list of GSEA objects (like what you get from gsea.split), extract
#' the geneset names, and create a gmx-like object. i.e. a table where each
#' column contains a vector of geneset names mentioned in each of the GSEA
#' objects in the list.
#' 
#' @param x a list of GSEA objects
#' @return a data.frame with colnames == names(x), no rownames, where each
#'   column contains the genesets which were in the relevant object inside x.
#' @author Mark Cowley, 2009-10-29
#' @export
gsealist2genesetlist <- function(x) {
	tmp <- lapply(x, function(x) names(x$leading.edge))
	N <- max(sapply(tmp, length))
	tmp <- lapply(tmp, function(x) c(x, rep(NA, N-length(x))))
	res <- as.data.frame(tmp)

	res
}
