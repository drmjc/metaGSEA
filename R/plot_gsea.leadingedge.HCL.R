#' GSEA Leading Edge HCL plot.
#' 
#' This function compares genesets in terms of the overlap of the leading edge
#' genes
#' within each geneset. This information is presented as a hierarchical
#' clustering plot,
#' where distance is measured proportional to the number of genes in common
#' (see \code{\link{jacquard}}).
#' The user can either pass in a leading edge list (ie a list where each
#' element is
#' a vector of gene symbols), or (as will be more common),
#' an entire leading edge list (eg 1434 elements for c2_all), in conjunction
#' with a GSEA
#' top table (tt), coupled with a threshold, such as the top 50, or those with
#' FDR < 0.25
#' via the N, FDR, P, or FWER arguments.
#' 
#' @note 
#'  max.distance\cr
#' The intuitive outcome of setting \code{max.distance} (eg 0.1)
#'   is that no genesets with a distance of >0.1 (ie a similarity of <90%)
#'   should remain in the HCL. This is not quite the case. We filter out genesets
#'   that have distance to all other genesets > 0.9 in the 'distance' object
#'   (see \code{\link{gsea.leadingedge.distance}}) which is slightly different to filtering
#'   out nodes in the resulting HCL (which is made upon the distance object).
#'   Net effect:\cr
#' - You can get singletons in the filtered HCL because all of
#'   the things that the geneset was connected to have been filtered out.\cr
#' - You can get genesets in a tree of significant genesets, which is itself
#'   below the distance threshold due to a similar explanation to above.
#' 
#' @param x a GSEA object. see \code{\link{import.gsea}}
#' @param N if not \code{NULL}, then choose the top N most significant genesets
#' @param NES if not \code{NULL}, then choose those genesets with |NES| > the
#'   threshold
#' @param FDR if not \code{NULL}, then choose those genesets with FDR < the
#'   threshold
#' @param P if not \code{NULL}, then choose those genesets with nominal P <
#'   the threshold
#' @param FWER if not \code{NULL}, then choose those genesets with FWER < the
#'   threshold
#' @param direction include just the genesets that are \dQuote{up}, \dQuote{down}, or
#'   \dQuote{either} means top 50 genesets that are either up or down
#' @param main additional arguments passed to plot.
#' @param xlab additional arguments passed to plot.
#' @param yaxis either \dQuote{similarity} or \dQuote{distance}. 
#'    Two closely related genesets have eg 90% similarity, and a distance of 0.1.
#' @param max.distance (experimental) attempt to filter out genesets that are
#'   dissimilar to all other genesets. This is a distance threshold, not a
#'   similarity threshold, so values of 0.9 are a good place to start. See
#'   note.
#' @param rename.genesets should the genesets be renamed to include the
#'   rank/fdr/direction
#' @param label.clusters if TRUE, then draw red boxes around the clusters. see
#'   cluster.threshold.
#' @param cluster.threshold the distance threshold if label.clusters=TRUE.
#'   0.99 is a good default in practice.
#' @param \dots additional arguments passed to plot.
#' @param hclust.method the hierarchical clustering method. see \code{\link{hclust}}
#' @return invisibly returns the hclust object.
#' @author Mark Cowley, 2009-04-06
#' @seealso \code{\link{plot_gsea.leadingedge.HCL}} \code{\link{plclust.gsea}} \code{\link{gsea.leadingedge.distance}}
#' @export
plot_gsea.leadingedge.HCL <- function(x,
	N=NULL, NES=NULL, FDR=NULL, P=NULL, FWER=NULL, direction=c("either", "up", "down"),
	main="Leading edge similarities", xlab="", 
	yaxis=c("similarity", "distance"),
	max.distance=1.0, rename.genesets=FALSE, 
	label.clusters=FALSE, cluster.threshold=0.99, 
	hclust.method="complete", ...) {

	stopifnot(is.gsea(x))
	
	direction <- direction[1]
	yaxis <- yaxis[1]
	MAX.LABEL.LENGTH <- 40

	#
	# subset the genesets from the leading edge by looking in the top table.
	#
	x <- gsea.filter(x, N=N, NES=NES, FDR=FDR, P=P, FWER=FWER, direction=direction)
	genesets <- x$tt$NAME

	#################
	# determine a suitable sub-title 
	# (do this before checking if there's enough genesets)
	sub <- ""
	if( !is.null(N) ) {
		N <- min(N, nrow(x$tt))
		sub <- paste("top", N)
	}
	else if( !is.null(NES) ) {
		sub <- paste("|NES| >", prettyNum(abs(NES)), paste(sep="", "(N=", length(genesets), ")"))
	}
	else if( !is.null(FDR) ) {
		sub <- paste("FDR <", prettyNum(FDR), paste(sep="", "(N=", length(genesets), ")"))
	}
	else if( !is.null(P) ) {
		sub <- paste("nominal P <", prettyNum(P), paste(sep="", "(N=", length(genesets), ")"))
	}
	else if( !is.null(FWER) ) {
		sub <- paste("FWER <", prettyNum(FWER), paste(sep="", "(N=", length(genesets), ")"))
	}
	else {
		sub <- "no filter"
	}
	sub <- paste(sub, direction, sep=" :: ")
	###################

	#
	# make sure there's enough genesets to do an HCL on...
	#
	if( length(genesets) <= 2 ) {
		plot.blank(main=main, sub=sub, box=TRUE, message="**** too few genesets pass thresholds ****")
		return(FALSE)
	}

	if( rename.genesets ) {
		#
		# augment the geneset name to include its rank and FDR.
		#
		x <- gsea.rename.genesets(x, rank=TRUE, fdr=TRUE, direction=(direction=="either"), maxlen=MAX.LABEL.LENGTH)
		# names(x$leading.edge) <- gsea.rename.genesets.deprecated(names(x$leading.edge), x$tt, direction=="either", MAX.LABEL.LENGTH)
	}

	#
	# calculate pairwise differences, and optionally exclude genesets that are dissimilar to all others??
	#
	d <- gsea.leadingedge.distance(x$leading.edge)
	if( max.distance < 1.0 ) {
		d <- threshold.distance.matrix(d, max.distance)
	}
	hc <- hclust(d, method=hclust.method)

	plclust.gsea(hc, main=main, xlab=xlab, sub=sub, yaxis="similarity", ...)

	# did we threshold the distance matrix? add a red line if we did.
	if( max.distance < 1.0 ) {
		abline(h=max.distance, col="red", lwd=2)
	}
	
	if( label.clusters ) {
		rect.hclust(hc, h=cluster.threshold)
		rect.hclust.labels(hc, h=cluster.threshold)
	}

	invisible(hc)
}


#' Take a GSEA object and calculate an HCL between the leading edge genes.
#' Don't worry about the fancy geneset renaming that plot_gsea.leadingedge.HCL
#' gives you.
#' 
#' @param x a GSEA object
#' @param main the plot title
#' @param xlab the x label
#' @param hclust.method the hierarchical clustering method. see \code{\link{hclust}}
#' @param \dots further arguments passed to \code{\link{plclust.gsea}}
#' @return none. creates an HCL dendrogram plot
#' @author Mark Cowley, 2009-09-01
#' @seealso \code{\link{plot_gsea.leadingedge.HCL}} \code{\link{plclust.gsea}} \code{\link{gsea.leadingedge.distance}}
#' @export
plot_gsea.leadingedge.HCL.simple <- function(x,
	main="Leading edge similarities", xlab="", hclust.method="complete", ...) {
	stopifnot(is.gsea(x))
	d <- gsea.leadingedge.distance(x$leading.edge)
	hc <- hclust(d, method=hclust.method)
	
	N <- nrow(hc$merge) + 1
	plclust.gsea(hc, main=main, xlab=xlab, yaxis="similarity", sub=paste(N, "genesets"), ...)
}


#' Take a GSEA object, do an HCL between the leading edge genes, cut the HCL, 
#' then plot 1 HCL per page of a PDF.
#' 
#' @param x a GSEA object
#' @param main the plot title
#' @param xlab the x label
#' @param hclust.method the hierarchical clustering method. see \code{\link{hclust}}
#' @param h numeric scalar or vector with heights where the tree should be
#'   cut. (see cutree)
#' @param \dots see \code{\link{plot_gsea.leadingedge.HCL.simple}}
#' @author Mark Cowley, 2009-09-01
#' @export
plot_gsea.leadingedge.HCL.1clusterPerPage <- function(x,
	main="Leading edge similarities", xlab="", hclust.method="complete", h=0.99, ...) {
	stopifnot(is.gsea(x))
	d <- gsea.leadingedge.distance(x$leading.edge)
	hc <- hclust(d, method=hclust.method)
	ct <- cutree(hc, h=h)
	
	for(i in sort(unique(ct))) {
		if( sum(ct==i) <= 3 ) {
			# plot.blank(main=paste("too few genesets in cluster", i))
			next
		}
		else {
			tmp <- subset.gsea(x, ct==i)
			plot_gsea.leadingedge.HCL.simple(tmp, main=paste(main, "-", "Cluster", i), xlab=xlab, hclust.method=hclust.method, ...)
		}
	}
}
