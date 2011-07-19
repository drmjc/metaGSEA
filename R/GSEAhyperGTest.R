#' Perform a Hypergeometric test, upon a GMT file.
#' This is a convenience wrapper around \code{\link{GMThyperGTest}}.
#' 
#' @param genes a vector of gene symbols
#' @param gmt a gmt object. see \code{\link{import.gsea.gmt}}
#' @param adjust.method the multiple hypothesis correction value. see \code{\link{p.adjust}}
#' @return a \code{data.frame} containing these columns: 
#' \dQuote{GeneSet}, \dQuote{x}, \dQuote{m}, \dQuote{n}, \dQuote{k}, \dQuote{Pval}, \dQuote{FDR},
#' where:\cr
#' \dQuote{x} = the number of white balls drawn without replacement from an urn which contains both black and white balls.\cr
#' \dQuote{m} = the number of white balls in the urn.\cr
#' \dQuote{n} = the number of black balls in the urn.\cr
#' \dQuote{k} = the number of balls drawn from the urn.\cr
#' @seealso \code{\link[stats]{phyper}}
#' @author Mark Cowley, 2009-01-15
#' @export
GSEAhyperGTest <- function(genes, gmt, adjust.method="none") {
	GMThyperGTest(genes=genes, gmt=gmt, adjust.method=adjust.method)
}


#' Perform a Hypergeometric test, upon a GMT file.
#' Much like \code{GOhyperGTest} from \code{topGO}, this function computes 
#' the overlaps between a vector
#' of gene id's of interest, and a collection of gene sets, stored within a gmt file
#' 
#' The type of identifier used doesn't matter much, as long as the same identifiers are
#' used in genes and gmt, eg Gene Symbols.
#' 
#' @param genes a vector of gene symbols
#' @param gmt a gmt object. see \code{\link{import.gsea.gmt}}
#' @param adjust.method the multiple hypothesis correction value. see \code{\link{p.adjust}}.
#'   default=\dQuote{none}
#' @return
#' a \code{data.frame} containing these columns: 
#' \dQuote{GeneSet}, \dQuote{x}, \dQuote{m}, \dQuote{n}, \dQuote{k}, \dQuote{Pval}, \dQuote{FDR}, where:\cr
#' \dQuote{x} = the number of white balls drawn without replacement from an urn which contains both black and white balls.\cr
#' \dQuote{m} = the number of white balls in the urn.\cr
#' \dQuote{n} = the number of black balls in the urn.\cr
#' \dQuote{k} = the number of balls drawn from the urn.\cr
#' @author Mark Cowley, 2009-01-15
#' @seealso \code{\link[stats]{phyper}}
#' @export
GMThyperGTest <- function(genes, gmt, adjust.method="none") {
	warning("I think there's a -1 error in the phyper command. Mark to check PINA code.")
	stopifnot(adjust.method %in% p.adjust.methods)
	
	N <- length(setdiff(unique(unlist(gmt)), c(NA, "---", "")))
	genes <- unique(genes)
	
	res <- as.data.frame(matrix(NA, nrow=length(gmt), ncol=7), stringsAsFactors=FALSE)
	# rownames(res) <- names(gmt)
	colnames(res) <- c("GeneSet", "x", "m", "n", "k", "Pval", "FDR")
	res[,1] <- names(gmt)
	
	for(gmt.idx in 1:length(gmt)) {
		# x	 the number of white balls drawn without replacement from an urn which contains both black and white balls.
		# m	 the number of white balls in the urn.
		# n	 the number of black balls in the urn.
		# k	 the number of balls drawn from the urn.
		x <- length(intersect(genes, gmt[[gmt.idx]]))
		m <- length(gmt[[gmt.idx]])
		n <- N - m
		k <- length(genes)
		hyper.result <- dhyper(x,m,n,k)
		res[gmt.idx, 2:5] <- c(x,m,n,k)
		res[gmt.idx, 6] <- hyper.result
	}
	
	res$FDR <- p.adjust(res$Pval, adjust.method)
	res <- res[order(res$Pval, decreasing=FALSE), ]
	rownames(res) <- 1:nrow(res)
	
	res
}
