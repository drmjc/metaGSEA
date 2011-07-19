# Take a GSEA list, and filter out genesets that fail various criteria.
# The criteria are additive.
#
# Parameters:
#	x: a gsea list
#	N: keep the top 'N' genesets, eg 50
#	NES: keep the genesets with |NES| > |NES threshold|.
#	FDR, P, FWER: keep the genesets with NOM.p.val, FDR.q.val or FWER.p.val < threshold
#	direction: "up", "down", or NULL or "either"
#	genesets: a vector of geneset ID's such that only these genesets will be retained.
#
# Value:
#	a GSEA list with probably fewer genesets in the top table, the leading edge list, and the edb object [if present].
#	also, a "filter.settings" list is added which details which settings were used
#
# Mark Cowley, 2009-09-03
# 2010-11-11: remove the =NULL and is.null check with calls to missing()
# 2010-11-16: reverse the =NULL change


##' Take a GSEA list, and filter out genesets that fail various criteria.
##' 
##' The criteria are additive.
##' 
##' @param x a gsea list
##' @param N keep the top 'N' genesets, eg 50
##' @param NES keep the genesets with |NES| > |NES threshold|.
##' @param FDR keep the genesets with NOM.p.val, FDR.q.val or FWER.p.val <
##'   threshold
##' @param P keep the genesets with NOM.p.val, FDR.q.val or FWER.p.val <
##'   threshold
##' @param FWER keep the genesets with NOM.p.val, FDR.q.val or FWER.p.val <
##'   threshold
##' @param direction "up", "down", or NULL or "either"
##' @param genesets a vector of geneset ID's such that only these genesets will
##'   be retained.
##' @return a GSEA list with probably fewer genesets in the top table, the
##'   leading edge list, and the edb object [if present]. also, a
##'   "filter.settings" list is added which details which settings were used
##' @author Mark Cowley, 2009-09-03
##' @export
gsea.filter <- function(x, N=NULL, NES=NULL, FDR=NULL, P=NULL, FWER=NULL, direction=NULL, genesets=NULL) {
	stopifnot(is.gsea(x))
	
	if( !missing(direction) && direction != "either" ) {
		if( direction == "up" )
			x$tt <- x$tt[x$tt$DIRECTION == "up", ]
		else if( direction == "down" ) {
			x$tt <- x$tt[x$tt$DIRECTION == "down", ]
			x$tt <- x$tt[nrow(x$tt):1, ]
		}
		else {
			stop("direction must be NULL, up or down")
		}
	}
	if( !is.null(N) ) {
		N <- min(N, nrow(x$tt))
		x$tt <- x$tt[order(abs(x$tt$NES), decreasing=TRUE), ]
		x$tt <- x$tt[1:N, ]
	}
	if( !is.null(NES) ) {
		x$tt <- x$tt[order(abs(x$tt$NES), decreasing=TRUE), ]
		x$tt <- x$tt[abs(x$tt$NES) > abs(NES), ]
	}
	if( !is.null(FDR) ) {
		idx <- which(x$tt$FDR.q.val < FDR)
		# x$tt <- x$tt[order(x$tt$FDR.q.val, decreasing=FALSE), ]
		x$tt <- x$tt[idx, ]
	}
	if( !is.null(P) ) {
		x$tt <- x$tt[order(x$tt$NOM.p.val, decreasing=FALSE), ]
		x$tt <- x$tt[x$tt$NOM.p.val < P, ]
	}
	if( !is.null(FWER) ) {
		x$tt <- x$tt[order(x$tt$FWER.p.val, decreasing=FALSE), ]
		x$tt <- x$tt[x$tt$FWER.p.val < FWER, ]
	}
	if( !is.null(genesets) ) {
		genesets <- intersect(genesets, x$tt$NAME)
		idx <- match(genesets, x$tt$NAME)
		x$tt <- x$tt[idx,]
	}
	
	x$leading.edge <- x$leading.edge[match(x$tt$NAME, names(x$leading.edge))]

	if( "edb" %in% names(x) && !is.na(x$edb) ) {
		x <- gsea.filter.edb(x, x$tt$NAME)
		x$gmt <- x$gmt[match(x$tt$NAME, names(x$gmt))]
	}
	x$filter.settings <- list(N=N, NES=NES, FDR=FDR, P=P, FWER=FWER, direction=direction, genesets=genesets)

	return( x )
}

#' Subsetting GSEA objects.
#' Return subsets of GSEA objects which meet conditions, with consistent args to
#' \code{\link[base]{subset}}
#' 
#' @param x object to be subsetted.
#' @param subset logical expression indicating elements or rows to keep: 
#'     missing values are taken as false.
#' @param \dots unused
#' @author Mark Cowley, 2011-02-21
#' @export
subset.gsea <- function(x, subset, ...) {
	if(!is.gsea(x)) stop("x must be a GSEA object.\n")
	if(!is.logical(subset)) stop("subset must be a logical vector.\n")
	if( length(subset) != length(x$leading.edge) ) stop("subset must be the same length as the number of genesets in x.\n")
	if( all(subset) ) return(x)
	ids <- names(x$leading.edge)[subset]
	res <- gsea.filter(x, genesets=ids)

	return(res)
}

# Filter a GSEA edb file, restricting the entries to restricted set of genesets.
#
# Parameters:
#	x: a gsea list (not just the edb file)
#	genesets: a vector of geneset names. Minimal error checking done, so make sure they exist in your edb file!
#
# Value:
#	a gsea list with the edb object restricted to only those that were specified by \code{genesets}
#
# Mark Cowley, 2009-10-13
# 2009-10-16: speed improvement by only iterating through relevant genesets.
#


##' Filter a GSEA edb file, restricting the entries to restricted set of
##' genesets.
##' 
##' @param x a gsea list (not just the edb file)
##' @param genesets a vector of geneset names. Minimal error checking done, so
##'   make sure they exist in your edb file!
##' @return a gsea list with the edb object restricted to only those that were
##'   specified by \code{genesets}
##' @author Mark Cowley, 2009-10-13
##' @export
gsea.filter.edb <- function(x, genesets) {
	stopifnot( "edb" %in% names(x) )
	
	require(XML)
	root <- xmlRoot(x$edb, skip=TRUE)

	edb.names <- .edb.names(x$edb)
	idx <- match(genesets, edb.names)

	suppressWarnings(res <- xmlTree("EDB", attrs=xmlAttrs( root )))
	for(i in idx) {
		a <- xmlAttrs(root[[i]])
		geneset <- sub("^.*#", "", a[3])
		if( geneset %in% genesets )
			res$addNode("DTG", attrs=a)
	}
	
	x$edb <- res
	
	return( x )
}
