#' Split a GSEA object based on leading edge similarities.
#' Take a GSEA object upon G genesets, and split it into k GSEA objects,
#' where all genesets within each object has similarity to all other genesets.
#' 
#' Internally, this function first runs \code{\link{gsea.leadingedge.distance}}, then
#' \code{\link{hclust}}, and finally \code{\link{cutree}}. The resulting cluster assignments
#' are then used to subet GSEA objects into smaller GSEA objects with only
#' specified genesets within them. You can then export them, and load each
#' subset in the LeadingEdge Viewer tool.
#' 
#' @param x a GSEA object
#' @param dist the distance to cut at. 0.99 is a large distance (close to the
#'   root of the tree).
#' @param min.size discard geneset clusters that have fewer than min.size
#'   genesets in them.
#' @return a named list of GSEA objects, one per cluster of genesets. The
#'   names are cluster001, cluster002, ...
#' @author Mark Cowley, 2009-10-16
#' @export
gsea.split <- function(x, dist=0.99, min.size=1) {
	d <- gsea.leadingedge.distance(x$leading.edge)
	hc <- hclust(d)
	tmp <- cutree(hc, h=dist)
	N <- max(tmp)
	cat(sprintf("Found %d clusters.\n", N))

	geneset.list <- list()
	for(i in 1:N) {
		geneset.list[[i]] <- names(tmp)[tmp==i]
	}
	names(geneset.list) <- sprintf("cluster%03d", 1:N)
	
	geneset.list <- geneset.list[sapply(geneset.list, length) >= min.size]
	cat(sprintf("%d clusters have min.size >= %d.\n", length(geneset.list), min.size))
	
	cat("Creating GSEA subsets")
	res <- list()
	for(i in 1:length(geneset.list)) {
		res[[i]] <- gsea.filter(x, genesets=geneset.list[[i]])
		cat(".")
	}
	cat("\n")
	names(res) <- names(geneset.list)

	return( res )
}
