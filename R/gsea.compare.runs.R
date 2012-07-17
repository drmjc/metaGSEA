#' Compare N GSEA runs vs a single gmt file (eg lots of rnk files vs c2_all)
#' 
#' @param dirs a vector of gsea directories
#' @param gsealist optional list of gsea runs that you've already imported.
#'   over-rides <dirs>
#' @param outfile if specified, an XLS output file will be created.
#' @param method which genesets to keep if they are not the same? \dQuote{intersect}:
#'   just the genesets found in all GSEA runs; \dQuote{union}: all the genesets found
#'   in any GSEA run. \code{NA}'s will be used for those genesets that were not
#'   observed.
#' @return a data.frame of \dQuote{geneset}, \dQuote{size}, \dQuote{NES}, \dQuote{FDR}, \dQuote{NES}, \dQuote{FDR}, \dots
#' @author Mark Cowley, 2009-12-16
#' @export
#' @importFrom excelIO write.xls
gsea.compare.runs.1gmt <- function(dirs, gsealist, outfile, method=c("intersect", "union")[1]) {
	if( !missing(gsealist) ) {
		gsea <- gsealist
	}
	else if( !missing(dirs) ) {
		gsea <- list()
		for(i in 1:length(dirs))
			gsea[[i]] <- import.gsea(dirs[i])
		
		# Name the gsea runs by their rnk file, or failing that, the rpt_label
		if( "rnk" %in% names(gsea[[1]]$rpt) )
			tmp.names <- sub("\\.rnk", "", basename(sapply(gsea, function(x) x$rpt$rnk)))
		else
			tmp.names <- sapply(gsea, function(x) x$rpt$rpt_label)
		names(gsea) <- tmp.names
	}
	else {
		stop("Must supply either dirs, or gsealist to gsea.compare.runs.1gmt.\n")
	}

	NES <- gsea.merge.tt(lapply(gsea, function(x) x$tt), keep="NES", method=method)
	FDR <- gsea.merge.tt(lapply(gsea, function(x) x$tt), keep="FDR.q.val", method=method)

	colnames(NES)[3:ncol(NES)] <- paste("NES.", colnames(NES)[3:ncol(NES)], sep="")
	colnames(FDR)[3:ncol(FDR)] <- paste("FDR.", colnames(FDR)[3:ncol(FDR)], sep="")

	res <- collate.data.frame(NES, FDR)
	res <- res[,-c(2,4)]

	# reorder the rows by decreasing average NES
	o <- order(rowMeans(NES[,3:ncol(NES)], na.rm=TRUE), decreasing=TRUE)
	# o <- order(res[,3], decreasing=TRUE)
	res <- res[o, ]

	if( !missing(outfile) )
		write.xls(res, outfile, na="")

	invisible(res)
}
