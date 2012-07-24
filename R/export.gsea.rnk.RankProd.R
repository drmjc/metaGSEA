#' Take the results from a RankProd analysis, and produce a rnk file for GSEA.
#' 
#' The score is the inverse of the RP score, multiplied by 100 or -100 for the
#' up and down-regulated genes respectively. This typically produces scores in
#' 3 to -3 range
#' 
#' @param rp a rankprod result object.
#' @param names a character vector of probe names, or if \code{NULL}, the rownames 
#'   of x will be used, else error
#' @param file the output rnk file name
#' @param method one of \sQuote{inv.log}, or \sQuote{inverse}. see \code{\link{RankProd.signed.score}}
#' @return write out a rnk file, and invisibly return the rp score.
#' @author Mark Cowley, 2009-01-09
#' @export
export.gsea.rnk.RankProd <- function(rp, names=NULL, file, method=c("inv.log", "inverse")) {
	require(microarrays) || stop("required package 'microarrays' is not installed")
	rp <- rp$RPs
	if( is.null(names) ) {
		if( is.null(rownames(rp)) )
			stop("Need names, or the rownames of the RPs to be specified.")
		else
			names <- rownames(rp)
	}
	score <- RankProd.signed.score(rp, method=method)
	
	export.gsea.rnk(score, names, file, TRUE)
	
	invisible(score)
}
