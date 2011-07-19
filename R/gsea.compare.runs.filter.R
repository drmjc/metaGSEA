#' Function to filter a set of GSEA run comparisons, and/or to re-order the
#' rows.
#' 
#' This also provides extensive ability to re-sort the filtered rows, via \code{sort.by} 
#'   and the \code{\link{gsea.compare.runs.sort}}?
#' 
#' If you want to just sort the rows, see \code{\link{gsea.compare.runs.sort}}, or set \code{fdr=NULL}.
#' 
#' @param x a data.frame from \code{\link{import.gsea.compare.runs}}
#' @param fdr a single numeric value, such as 0.05 or NULL
#' @param mode should \dQuote{any}, or \dQuote{all} of the values satisfy the filter?
#' @param sort.by one of \dQuote{none}, \dQuote{avgabs}, \dQuote{range}, \dQuote{deltaAB}, or an integer. see \code{\link{gsea.compare.runs.sort}}
#' @return x with only the rows that satisfy the filter.
#' @author Mark Cowley, 2009-03-23
#' @export
gsea.compare.runs.filter <- function(x, fdr=NULL, mode=c("any", "all")[1], sort.by=c("none", "avgabs", "range", "deltaAB", 1)) {
	sort.by <- sort.by[1]
	
	#
	# filter out rows that fail criteria
	#
	idx <- 1:nrow(x)
	if( !is.null(fdr) ) {
		FDR <- x[, grep("FDR", colnames(x))]
		colnames(FDR) <- sub("FDR.q.val.", "", colnames(FDR))
		if( mode == "any" ) {
			rmin <- apply(FDR, 1, min, na.rm=TRUE)
			idx <- which(rmin < fdr)
		}
		else if( mode == "all" ) {
			nc <- ncol(FDR)
			rsum <- rowSums(FDR < fdr)
			idx <- which(rsum == nc)
		}
		else {
			stop("mode must be one of 'any' or 'all'")
		}
	}
	x <- x[idx, ]
	#
	# at this point, x has been filtered.
	#
	
	x <- gsea.compare.runs.sort(x, sort.by)

	return( x )
}
