#' Plot an hclust object, that has been made from GSEA genesets.
#' 
#' @param hc an object of class \code{hclust}
#' @param main see \code{\link{par}}
#' @param xlab see \code{\link{par}}
#' @param sub see \code{\link{par}}
#' @param yaxis what type of yaxis scale do you want? \dQuote{similarity}, or \dQuote{distance}?
#' @param lab.col argument to change the colour of the text labels. (not yet
#'   implemented)
#' @param hang see \code{\link{plclust}} default=0.03
#' @param cex see \code{\link{par}}. default=NULL
#' @param cex.main see \code{\link{par}}. default=1
#' @param \dots arguments passed to \code{\link{plclust}}'
#' @return none. creates a plot
#' @seealso \code{\link{plclust}}
#' @author Mark Cowley, 2009-09-04
#' @export
#' @rdname plclust.gsea
plclust.gsea <- function(hc, main="Leading edge similarities", xlab="", sub=NULL, yaxis=c("similarity", "distance"), lab.col="black", hang=0.03, cex=NULL, cex.main=1, ...) {
	
	yaxis <- yaxis[1]
	
	N <- length(hc$height) + 1
	
	if( is.null(cex) ) {
		if( N < 100 )
			cex <- 0.5
		else if( N < 150 )
			cex <- 0.35
		else 
			cex <- 0.25
	}

	lab.col <- recycle(lab.col, length(hc$order))

	plot(hc, main=main, ylab="", xlab=xlab, frame.plot=TRUE, hang=hang, axes=FALSE, cex=cex, sub=sub, cex.main=cex.main, ...)
	# plclust_in_colour(hc, main=main, ylab="", xlab=xlab, frame.plot=TRUE, hang=hang, axes=FALSE, cex=cex, sub=sub, lab.col=lab.col, ...)
	# plclust_in_colour2(hc, main=main, ylab="", xlab=xlab, frame.plot=TRUE, hang=hang, axes=FALSE, cex=cex, sub=sub, col=lab.col, ...)
	abline(h=seq(0,0.5,0.1), lty="dashed", col="grey")
		
	#
	# add an axis
	#
	if( yaxis == "similarity" ) {
		axis(side=2, at=seq(0,1,0.1), seq(100,0,-10))
		mtext(side=2, outer=TRUE, "Similarity (%)", line=-1, las=0)
	}
	else if( yaxis == "distance" ){
		axis(side=2, at=seq(0,1,0.1), seq(0,1,0.1))
		mtext(side=2, outer=TRUE, "Distance (1-Jacquard)", line=-1, las=0)
	}
}

#' @export
#' @rdname plclust.gsea
plot_gsea.hclust <- plclust.gsea
