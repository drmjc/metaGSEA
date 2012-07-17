#' Merge two GSEA objects.
#' 
#' You would want to do this if you compare two rnk files to 1 gmt file (eg
#' c2_all), and want to compare the significant genesets between these runs.
#' It's no good knowing that geneset X is in the top 5 in both runs if they
#' contain markedly different genes in their leading edges.
#' 
#' @param gsea1 GSEA object to be merged.
#' @param gsea2 GSEA object to be merged.
#' @param prefix character(2) containing optional prefix; or \code{NULL}
#' @param suffix character(2) containing optional suffix; or \code{NULL}
#' @param warn logical: produce a warning if \code{prefix=NULL && suffix=NULL}
#' @author Mark Cowley, 2009-10-14
#' @export
gsea.merge <- function(gsea1, gsea2, prefix=NULL, suffix=NULL, warn=TRUE) {
	if( gsea1$rpt$chip != gsea2$rpt$chip ) {
		# different chips, or different version of the same chip will invariably have different numbers of genes in them...
		# ... leading to different numbers of rows in the collapsed rnk files...
		# ... which will kill the Leading Edge viewer for the down-regulated sets
		# (since these sets contain negative genes, and the most negative genes tend 
		#  to be at the bottom of the collapsed rnk file.)
		stop( sprintf("You have to use the same .chip file to generate your GSEA results. You've used: %s and %s\n", shQuote(gsea1$rpt$chip), shQuote(gsea2$rpt$chip)) )
	}
	if( length(gsea1$rnk) != length(gsea2$rnk) ) {
		# for the same reason as above, different lengthed rnk files will kill the leading edge viewer tool.
		stop("You have different lengthed rnk files. You should run GSEA using all probesets on the chip, using the same .chip file, so that you get .rnk files that are the same length.\n")
	}


	# rename the gsea genesets.
	if( is.null(prefix) && is.null(suffix) ) {
		if( warn ) warning("I hope you have already renamed your gsea runs with gsea.rename.genesets! Otherwise you can't tell the origin of a geneset if it's listed twice.")
	}
	else {
		gsea1 <- gsea.rename.genesets(gsea1, prefix=prefix[1], suffix=suffix[1])
		gsea2 <- gsea.rename.genesets(gsea2, prefix=prefix[2], suffix=suffix[2])
	}

	# merge the tt's, and re-sort by NES
	tt <- rbind(gsea1$tt, gsea2$tt)
	tt <- tt[order(tt$NES, decreasing=TRUE), ]
	
	# merge the leading.edge's
	le <- lbind(gsea1$leading.edge, gsea2$leading.edge)
	le <- le[tt$NAME]
	
	res <- list(tt=tt, leading.edge=le, rnk1=gsea1$rnk, rnk2=gsea2$rnk, rpt1=gsea1$rpt, rpt2=gsea2$rpt)

	# merge the edb's
	if( "edb" %in% names(gsea1) ) {
		res$edb <- gsea.merge.edb(gsea1$edb, gsea2$edb)
		res$gmt <- lbind(gsea1$gmt, gsea2$gmt)
		res$gmt <- res$gmt[tt$NAME]
	}
	
	return( res )
}



#' Merge two GSEA edb objects.
#' 
#' @param edb1 two edb XML objects. They should already have been renamed to
#'   make the genesets unique. See gsea.rename.genesets
#' @param edb2 two edb XML objects. They should already have been renamed to
#'   make the genesets unique. See gsea.rename.genesets
#' @return an edb XML object with as many elements as the two input edb's
#'   combined.
#' @author Mark Cowley, 2009-10-14
#' @export
gsea.merge.edb <- function(edb1, edb2) {
	r1 <- xmlRoot(edb1, skip=TRUE)
	r2 <- xmlRoot(edb2, skip=TRUE)
	
	a <- NULL
	suppressWarnings(res <- xmlTree("EDB", attrs=xmlAttrs( r1 )))
	for(i in 1:xmlSize(r1)) {
		a <- xmlAttrs(r1[[i]])
		res$addNode("DTG", attrs=a)
	}
	rnk <- a[1]
	for(i in 1:xmlSize(r2)) {
		a <- xmlAttrs(r2[[i]])
		a[1] <- rnk # GSEA complains if you have different rnk files present.
		res$addNode("DTG", attrs=a)
	}
	
	res
}
