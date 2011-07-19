#' Apply Fisher's method to a GSEA list
#' 
#' Fisher's method is a form of meta-analysis, where a single p-value is
#' computed for a set of input p-values. This uses the p-values within the
#' top table objects.
#' The GSEA results must have been run against the same GMT file.
#' 
#' @param gsea.list a list of GSEA objects
#' @param values Which statistic to use, the unadjusted P-values (\dQuote{NOM.p.val}), 
#'    or the FDR q values (\dQuote{FDR.q.val})?
#' @param direction either: ignore the direction of the gene's change.
#' @return a data.frame with ID, {S, p.value, q.value, rank} all produced by
#'   fisherMethod, and the p-values from each of the N GSEA results.
#' @author Mark Cowley, 2010-07-22
#' @export
fisherMethod.GSEAlist <- function(gsea.list, values=c("NOM.p.val", "FDR.q.val"), direction=c("either", "same")) {
	require(MADAM)
	
	stopifnot(is.gsea.list(gsea.list))
	
	#
	# check that all results were compared to the same gmt/gmx file, which is
	# reported in the rpt
	#
	gmx <- n <- sapply(gsea.list, function(x) x$rpt$gmx)
	if( !alleq(gmx) ) {
		stop(paste("The GseaPreRanked results must all have been generated using the same gmx/gmt file (eg c2.all.v2.5.symbols.gmt). We found these classes:", gmx, sep="\n"))
	}

	# make sure the specified 'values' are a valid column within the GSEA toptable.
	values <- values[1]
	stopifnot(values %in% colnames(gsea.list[[1]]$tt))

	pvals <- gsea.merge.tt(gsea.list, keep=values)
	colnames(pvals)[3:ncol(pvals)] <- paste(colnames(pvals)[3:ncol(pvals)],values, sep=".")
	
	tmp <- fisherMethod(pvals[,3:ncol(pvals)])
	res <- cbind(pvals, tmp)
	res <- res[order(res$rank), ]
	rownames(res) <- res$rank
	res$rank <- NULL
	
	res
}
