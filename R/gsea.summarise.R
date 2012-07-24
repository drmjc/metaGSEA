#' Summarise GSEA results
#' 
#' Summarise a single GSEA result, or list of results in terms of how many
#' genesets are changed.
#' 
#' @param x either a GSEA object, or a list of GSEA objects.
#' @return a \code{vector} or \code{data.frame} of number of gene sets passing various
#'   statistical thresholds See also: \code{/path/to/metaGSEA/bin/gsea.compare.runs.plot.sh}
#' @author Mark Cowley, 2009-03-23
#' @export
gsea.summarise <- function(x) {
	if( is.list(x) && ! "tt" %in% names(x) ) {
		res <- lapply(x, gsea.summarise)
		res <- cbind.list(res)
		res
	}
	else {
		res <- c(
			sum(x$tt$DIRECTION == "up"),
			sum(x$tt$DIRECTION == "down"),
			sum(x$tt$FDR.q.val < 0.25 & x$tt$DIRECTION == "up"), 
			sum(x$tt$FDR.q.val < 0.25 & x$tt$DIRECTION == "down"), 
			sum(x$tt$NOM.p.val < 0.01 & x$tt$DIRECTION == "up"), 
			sum(x$tt$NOM.p.val < 0.01 & x$tt$DIRECTION == "down"), 
			sum(x$tt$NOM.p.val < 0.05 & x$tt$DIRECTION == "up"),
			sum(x$tt$NOM.p.val < 0.05 & x$tt$DIRECTION == "down") 
		)
		names(res) <- 
			paste(rep(c("NUMBER", "FDR < 25%", "P < 0.01", "P < 0.05"), each=2), c("UP", "DOWN"), sep=" ")
		res
	}
}
