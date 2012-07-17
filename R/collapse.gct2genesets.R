#' Collapse a GCT file from 1 row per gene to 1 row per geneset.
#' 
#' This function collapses gene data to geneset data.
#' NB: the gct[,1] must be the same as the ID's in the gmt file. Usually, GCT
#' files have probe names in column 1, and GMT
#' files have gene symbols. Thus you need to collapse the GCT file to 1 row
#' per gene symbol (see collapse.gct).
#' Currently supported methods include:
#' mean, median which are good if you expect all genes in the pathway to be
#' similarly changed
#' min, max if you're looking at pathways with high or low abundance
#' sum which is good if the gct file is an adjacency matrix (eg summing up the
#' number of genes that are mutated in a pathway)
#' 
#' @param gct.file the path to a gct file
#' @param gmt.file the path to a gmt file. untested with a gmx file
#' @param method one of mean, median, min, max, sum
#' @param min.size minimum size threshold for each geneset. default=1, ie
#'   include all genesets with >= 1 gene. GSEA default = 15
#' @param max.size maximum size threshold for each geneset. default=1e05, ie
#'   include all genesets. GSEA default = 500
#' @param gct.out the path to the output gct file. This will overwrite any
#'   previous file there.
#' @return side effect of creating a file at gct.out, with 1 row per geneset,
#'   as long as the geneset passed the size thresholds.  Todo: add PGSEA,
#'   maxmean, maxvar as an option See also: collapse.gct
#' @author Mark Cowley, 2011-02-27
#' @export
collapse.gct2geneset <- function(gct.file, gmt.file, method=c("mean", "median", "min", "max", "sum"), min.size=1, max.size=1e05, gct.out=NULL) {
	if( is.null(gct.out) ) stop("Must specify the gct.out file.\n")
	if( !file.exists(gct.file) ) stop("gct.file does not exist.\n")
	if( !file.exists(gmt.file) ) stop("gmt.file does not exist.\n")
	if( min.size < 1 ) { cat("min.size is < 1. Forcing it to 1.\n"); min.size <- 1}

	gct <- import.gsea.gct(gct.file)
	gmt <- import.gsea.gmt(gmt.file)

	all.genes <- sort(unique(unlist(gmt)))
	if( length(intersect(gct[,1], all.genes)) < 0.01*nrow(gct) ) {
		stop("Very few gene symbols in common b/w gct and gmt file. Make sure they both use the same type of ID's. You probably need to collapse your gct file from probes to gene symbols. see collapse.gct\n")
	}

	# filter to genes in the gct, and then on size.
	gmt <- lapply(gmt, function(x) intersect(x, gct[,1]))
	gmt.len <- sapply(gmt, length)
	if( any(gmt.len < min.size) | any(gmt.len > max.size) ) {
		idx <- (gmt.len >= min.size) & (gmt.len <= max.size)
		gmt <- subset(gmt, idx)
		msg <- sprintf("After resticting genesets (a) to those present in the gct file, and (b) to those with at least %d genes, and (c) at most %d genes, leaves %d genesets.\n", min.size, max.size, length(idx))
		cat(msg)
	}
	if( length(gmt)==0 ) {
		msg <- "All genesets filtered out.\n"
		cat(msg)
		return()
	}
	
	# collapse the genes to genesets via various methods.
	res <- matrix(0, length(gmt), ncol(gct)-2)
	rownames(res) <- names(gmt); colnames(res) <- colnames(gct)[3:ncol(gct)]
	if( method %in% c("sum", "mean", "median", "min", "max") ) {
		for(i in 1:length(gmt)) {
			res[i,] <- apply(gct[match(gmt[[i]], gct[,1]), 3:ncol(gct)], 2, get(method), na.rm=TRUE)
		}
	}
	else {
		stop("Unsupported collapse method.\n")
	}
	
	desc <- names(gmt); names(desc) <- names(gmt)
	export.gsea.gct(res, description=desc, file=gct.out)
}
