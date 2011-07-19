#' Function to sort the rows of a GSEA comparison object
#' 
#' This provides extensive ability to re-sort the rows, via \code{sort.by}.\cr
#' \dQuote{avgabs} means sort by decreasing avg absolute NES, so genesets with large NES's in all comparisons will be
#'   prioritised;\cr
#' \dQuote{range}: genesets with large swings in NES across \emph{all} GSEA runs will be prioritised\cr
#' \dQuote{deltaAB}: genesets are ranked by the magnitude of the difference
#'   between \emph{only} run 1 and 2; or\cr
#' supply a single integer to specify which of the N runs you want to sort by. 
#' 
#' @param x a data.frame from \code{\link{import.gsea.compare.runs}}
#' @param sort.by one of \dQuote{none}, \dQuote{avgabs}, \dQuote{range}, \dQuote{deltaAB}, or an integer. see details
#' @return x with the rows resorted
#' @author Mark Cowley, 2011-07-19
#' @seealso \code{\link{gsea.compare.runs.filter}}
#' @export
gsea.compare.runs.sort <- function(x, sort.by=c("none", "avgabs", "range", "deltaAB", 1)) {
	sort.by <- sort.by[1]

	#
	# re-sort the rows
	#
	new.order <- 1:nrow(x)

	NES <- x[, grep("NES", colnames(x))]
	colnames(NES) <- sub("NES.", "", colnames(NES))

	if( sort.by == "none") {
		# do nothing
	}
	else if( sort.by == "avgabs") {
		absNES <- abs(NES)
		avgNES <- rowMeans(absNES)
		new.order <- order(avgNES, decreasing=TRUE)
	}
	else if ( sort.by == "range" ) {
		rangeNES <- apply(NES, 1, function(x) diff(range(x)))
		new.order <- order(rangeNES, decreasing=TRUE)
		
	}
	else if( sort.by == "deltaAB" ) {
		delta <- NES[,1] - NES[,2]
		# delta <- rep(0, nrow(NES))
		# for(i in 1:length(delta)) {
		# 	delta[i] <- NES[i,1] - NES[i,2]
		# }
		new.order <- order(delta, decreasing=TRUE)
	}
	else if ( is.numeric(sort.by) ) {
		new.order <- order(NES[,sort.by], decreasing=TRUE)
	}
	
	x <- x[new.order, ]
	
	return( x )
}
