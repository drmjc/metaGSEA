#' Convert a GSEA list into an adjacency matrix at the gene-level.
#' 
#' From a list of GSEA objects (like what you get from gsea.split), extract
#' the genes from the leading edges of each GSEA object, and then convert into
#' an adjacency matrix.
#' 
#' @param x a list of GSEA objects.
#' @return a \code{data.frame} of 0 or 1 depending on whether each gene was found in
#'   each element of the GSEA list.
#' 
#' @seealso \code{\link[mjcbase]{list2adjmat}} \code{\link{gsea2adjmat}}
#' @author Mark Cowley, 2009-10-29
#' @export
#' 
gsealist2adjmat <- function(x) {
	x <- sapply(x, "[", "leading.edge")
	names(x) <- sub(".leading.edge", "", names(x))
	x <- lapply(x, unlist)
	
	res <- list2adjmat(x)
	return(res)
}


#' Convert a GSEA object into an adjacency matrix at the gene-level.
#' 
#' From a GSEA object, extract the genes from the leading edges of the
#' genesets, and then convert into an adjacency matrix.
#' 
#' @param x a GSEA list
#' @return a data.frame of 0 or 1 depending on whether each gene was found in
#'   each element of the GSEA list.  
#' 
#' @seealso \code{\link[mjcbase]{list2adjmat}} \code{\link{gsealist2adjmat}}
#' @author Mark Cowley, 2009-10-29
#' @export
gsea2adjmat <- function(x) {
	res <- list2adjmat(x$leading.edge)

	return(res)
}
