#' Merge N GSEA toptables, retaining a column of interest (eg NOM.p.val, or
#' FDR...)
#' 
#' @param \dots a single gsea.list object, a list of gsea toptables, or
#'   individually named GSEA toptables.
#' @param keep which column to keep: \dQuote{NAME}, \dQuote{SIZE}, \dQuote{RANK}, \dQuote{DIRECTION},
#'   \dQuote{ENRICHED.PHENOTYPE}, \dQuote{ES}, \dQuote{NES}, \dQuote{NOM.p.val}, \dQuote{FDR.q.val}, \dQuote{FWER.p.val},
#'   \dQuote{RANK.AT.MAX}, \dQuote{LEADING.EDGE}, \dQuote{LEADING.EDGE.SIZE}
#' @param method which genesets to keep if there are differences across runs
#'   \dQuote{intersect} - just the genesets found in all GSEA runs; \dQuote{union} - all the
#'   genesets found in any GSEA run. \code{NA}'s will be used for those genesets that
#'   were not observed.
#' @author Mark Cowley
#' @export
gsea.merge.tt <- function(..., keep="NES", method=c("intersect", "union")[1]) {
	if( missing(keep) )
		stop("Need to specify which column to keep, and it must be named. eg keep=\"NES\"\n")
	else
		keep <- keep[1]
	
	tt.list <- list(...)
	if( length(tt.list) == 1 && is.gsea.list(tt.list[[1]]) ) {
		n <- names(tt.list[[1]])
		tt.list <- lapply(tt.list[[1]], function(x) x$tt)
		names(tt.list) <- n
	}
	else if( length(tt.list) == 1 && is.list(tt.list[[1]]) && !is.data.frame(tt.list[[1]]) ) {
		tt.list <- tt.list[[1]]
	}

	geneset.names <- sapply(tt.list, "[", "NAME")
	
	# genesets <- sort(unique(unlist(geneset.names)))
	if(method=="intersect") {
		genesets <- do.call("intersectN", args=geneset.names)
	}
	else if(method=="union") {
		# genesets <- sort(union.list(geneset.names))
		genesets <- do.call("unionN", args=geneset.names)
	}
	genesets <- sort(genesets)
	
	pvals <- matrix(NA, nrow=length(genesets), ncol=length(tt.list))
	colnames(pvals) <- names(tt.list)
	#
	# use 1 matrix to remember the 'keep' values, and 1 to remember the geneset SIZE
	#
	sizes <- res <- data.frame(NAME=genesets, SIZE=rep(0, length(genesets)), pvals)
	for(i in 1:length(tt.list)) {
		tt.list[[i]] <- tt.list[[i]][match(genesets, tt.list[[i]]$NAME), ]
		# idx <- match(tt.list[[i]]$NAME, genesets)
		res[,i+2] <- tt.list[[i]][,match(keep, colnames(tt.list[[i]]))]
		sizes[,i+2] <- tt.list[[i]][,match("SIZE", colnames(tt.list[[i]]))]
	}
	res$SIZE <- pwbc::rowMax(sizes[,3:ncol(sizes)], na.rm=TRUE)
	
	res$NAME <- as.character(res$NAME)
	
	res
}
