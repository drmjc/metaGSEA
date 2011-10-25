#' plot a panel of gsea histograms
#' From a GSEA object, plot a histogram of the nominal p values,
#' FDR q values, or FWER p values.
#'
#' @param x a gsea object
#' @param main the plot titles
#' @param which one or more of \dQuote{p}, \dQuote{q}, or \dQuote{FWER}, 
#'  corresponding to these columns in a GSEA report: \dQuote{NOM.p.val}, 
#'  \dQuote{FDR.q.val}, \dQuote{FWER.p.val}, respectively.
#' @return none. creates a panel of histograms.
#' @author Mark Cowley, 2011-10-20
#' @export
plot.gsea.statistic.histograms <- function(x, main="", which=c("p", "q", "FWER")) {
	!missing(x) && is.gsea(x) || stop("x must be a gsea object")
	which <- match.arg(which, c("p", "q", "FWER"), several.ok=TRUE)
	
	if( length(which) > 1 ) par(mfrow=c(1,length(which)))

	col <- colours.mjc("reds", 20)
	for(i in 1:length(which)) {
		cn <- switch(which[i],
			p="NOM.p.val",
			q="FDR.q.val",
			FWER="FWER.p.val"
		)
		vals <- x$tt[,cn]
		hist(vals, breaks=seq(0,1,0.05), main=main, xlab=cn, col=col)
	}
}
