#' Export a res file for GenePattern.
#' 
#' This requires a matrix of data (such as from RMA) and a matrix of calls,
#' usually MAS5 calls
#' of P/M/A. If you use Affymetrix ST arrays, you can generate DABG calls,
#' then use \code{\link[mjcaffy]{dabg2calls}}
#' 
#' @param data matrix of expression levels. If this is log base 2, then it
#'   will be unlog2ged. If you use any other log base, then unlog2 it yourself
#' @param calls matrix of calls: P/M/A
#' @param description a description for each ProbeSet. If NULL, then rownames
#'   of data will be used. if it doesn't have names, we assume it is the same
#'   length and in samre order as nrow(data) if it does have names, then only
#'   the relevant values will be extracted using the rownames(data) as the
#'   key.
#' @param file the filename
#' @param missing the value to write if there's any NA values
#' @param unlog2 logical: unlog the log-base-2 \code{data}? 
#' @return creates a res file
#' @author Mark Cowley, 2009-06-26
#' @seealso \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#res}
#' @export
export.broad.res <- function(data, calls, description=NULL, file, missing="", unlog2=TRUE) {
	if( is.null(description) ) {
		description <- rep("", nrow(data))
		names(description) <- rownames(data) # this avoids the warning in next section
	}

	if( !is.null(names(description)) ) {
		description <- description[rownames(data)]
	}
	else {
		warning("We're assuming that the order of rows in data, and description are the same.\n")
	}
	
	if( unlog2 && max(data, na.rm=TRUE) < 32 ) {
		data <- 2^data
		data <- round(data, 1)
	}
	else {
		data <- round(data, 4)
	}

	res <- affy.pivot.table(data, calls)
	res <- as.data.frame(res, stringsAsFactors=FALSE)
	res <- rownames2col(res, 1, "Accession")		
	res$Description <- description
	res <- move.column(res, "Description", 1)
	
	line1 <- colnames(res)
	line1[grep("Detection", line1)] <- ""
	line1 <- sub("_Signal", "", line1)
	line1 <- line1[1:(length(line1)-1)]
	line1 <- paste(line1, collapse="\t")
	
	tmp <- matrix("\t", 3, ncol(data))
	tmp[2,] <- colnames(data)
	line2 <- as.vector(tmp) # byrow=FALSE
	line2 <- line2[1:(length(line2)-1)]
	line2 <- paste(line2, collapse="")
	
	line3 <- nrow(res)
	
	OUT <- file(file, "w")
	writeLines(c(line1,line2,line3), OUT)
	write.delim(res, OUT, row.names=FALSE, col.names=FALSE, na=missing)
	close(OUT)
}
# CHANGELOG:
# 2013-01-17: default description is now ""

#' merge expression and calls together
#' 
#' @param exp a \code{matrix}-like object of expression levels
#' @param calls a \code{matrix}-like object of detection calls (P/M/A)
#' @return none.
#' @author Mark Cowley
#' @export
affy.pivot.table <- function(exp, calls) {
	# if( !is.matrix.like(exp) ) exp <- exprs(exp)
	# if( !is.matrix.like(calls) ) calls <- exprs(calls)
	exp <- round(exp,4)
	
    res <- interleave.columns(exp, calls)
    colnames(res) <- sub(".CEL", "", colnames(res))
    tmp2 <- rep(c("_Signal", "_Detection"), ncol(exp))
    colnames(res) <- paste(colnames(res), tmp2, sep="")

    invisible(res)
}
# CHANGELOG:
# 2013-01-16: 
# - dropped calls to exprs
# - dropped file argument
