#' From a list of GSEA's, find the elements that contain particular genesets
#' of interest.
#' 
#' @section Scenario:
#' You've used \code{\link{gsea.split}} to split a LARGE GSEA run into ~100 smaller, related
#' GSEA runs. From a small list of genesets of interest (represented by
#' patterns), you want to know which of these smaller GSEA runs contain your
#' genesets of interest.
#' 
#' @param gsea.list a list of GSEA objects.
#' @param patterns a vector of patterns which should be unique enough to find
#'   just one cluster
#' 
#' @return a vector of geneset indices where each pattern was found. IF one of
#'   the patterns is not unique enough, then a list is returned.
#' 
#' @author Mark Cowley, 2009-10-16
#' @export
gsealist.find.geneset <- function(gsea.list, patterns) {
	gsea.names <- lapply(gsea.list, function(gsea) gsea$tt$NAME)
	res <- list()
	for(pattern in patterns) {
		res[[pattern]] <- lgrep(pattern, gsea.names, fixed=TRUE)
	}
	
	if( all(sapply(res, length)==1) ) {
		res <- unlist(res)
		names(res) <- names(gsea.list)[res]
	}
	else if( any(sapply(res, length) > 1) ) {
		warning("Your patterns were not specific enough! I am returning a list.")
	}
	else if( any(sapply(res, length) == 0) ) {
		warning("some patterns were not found!")
	}

	res
}
