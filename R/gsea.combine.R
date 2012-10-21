#' Combine N gsea runs together
#' 
#' This merges N independent GSEA results into a single much larger GSEA result
#' object. So that you don't end up with N copies of each geneset, you need to
#' rename the genesets before running this code to uniquify them, eg adding
#' a short prefix tag (eg \dQuote{MCF7} or \dQuote{T47D}) is enough to uniquify the
#' geneset names.
#' You should probably also run some sort of geneset filtering, to keep the result
#' size manageble, and since filtering on this combined 
#' list may fail.
#' 
#' @section TODO:
#'  delete the gsea.runs param; users can just use \code{call("gsea.combine", my.gsea.list)}
#'
#' @param \dots at least 2 individual GSEA objects, or
#' @param gsea.runs a list of GSEA objects. This overrides \dots
#' @param verbose logical.
#' @return a GSEA object with larger tt, and larger leading.edge objects.
#'   Currently, the rnk and edb objects are set to NA -- WARNING: this may
#'   make routines that work on GSEA objects fail. We recommend treating
#'   results from gsea.combine as immutable objects.
#' @author Mark Cowley, 2010-08-16
#' @seealso \code{\link{gsea.rename.genesets}} \code{\link{gsea.filter}}
#' @export
gsea.combine <- function(..., gsea.runs, verbose=TRUE) {
	if( missing(gsea.runs) )
		gsea.runs <- list(...)

	if( length(gsea.runs) < 2 ) {
		stop("You must combine more than one GSEA run together. If you already have a list of GSEA runs, specify the gsea.runs arg, otherwise, provide each individual GSEA object.\n")
	}
	
	res <- list()
	res <- gsea.runs[[1]]
	res$rnk <- NA
	res$edb <- NA
	if( verbose ) message("Method cannot combine res and edb objects - these will be discarded.\n")
	for(i in 2:length(gsea.runs)) {
		res$tt <- rbind(res$tt, gsea.runs[[i]]$tt)
		res$leading.edge <- lbind(res$leading.edge, gsea.runs[[i]]$leading.edge)
		res$rnk
	}
	
	res
}
# CHANGELOG
# 2012-10-20: changed 'as.list(...)' to 'list(...)'
