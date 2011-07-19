#' Import a GenePattern res file
#' 
#' res files contain 2 columns per sample: a Signal, and a Call column, where Call contains
#' a detection call - either a P/M/A, or a proabability. These calls can be imported, or
#' thrown away via the \code{calls} parameter.
#' 
#' @param f the path to a res file
#' @param calls logical: if \code{TRUE}, retain the detection calls, if \code{FALSE}, throw away Calls, 
#' and rename Signal columns to strip out the trailing \dQuote{.Signal}
#' @return a \code{data.frame} of usually gene expression data. columns 1&2 are
#'   annotation, 3+ are data, in alternating pairs of Signal and Calls... ie
#'   2n+2 columns in total. If \code{calls=TRUE}, ncol=2+N, if \code{calls=FALSE}, ncol=2+2N
#' @author Mark Cowley, 2009-07-27
#' @export
import.gsea.res <- function(f, calls=TRUE) {
	header <- trim(readLines(f,3))
	nrow <- as.numeric(header[3])

	res <- read.delim(f, skip=3, header=FALSE, stringsAsFactors=FALSE)
	if( nrow != nrow(res) )
		stop(sprintf("You specified %d genes, but we found %d genes.\n", nrow, nrow(res)))

	tmp <- strsplit(header[1], "\t")[[1]]
	sample.names <- tmp[seq(3,ncol(res), 2)]
	colnames <- c(c("Description", "Accession"), paste(rep(sample.names, each=2), c("Signal", "Call"), sep="."))
	colnames(res) <- colnames
	rownames(res) <- res$Accession

	if( !calls ) {
		res <- res[,-grep("Call$", colnames(res))]
		colnames(res) <- sub("\\.Signal", "", colnames(res))
	}

	res
}
