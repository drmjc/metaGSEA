#' test whether a connection to BROAD's FTP/HTTP sites can be
#' made.
#' 
#' @param gsea.url a URL to test
#' @return logical(1) with an error message printed to stdout if no connection
#'   can be made.
#' @author Mark Cowley, 2009-12-10
#' @export
test.broad.ftp <- function(gsea.url) {
	con <- url(gsea.url)
	tmp <- NULL
	try(tmp <- readLines(con, 1), silent=TRUE)
	close(con)
	
	success <- FALSE
	if( !is.null(tmp) ) {
		success <- TRUE
	}
	else {
		url.that.should.exist <- "ftp://gseaftp.broad.mit.edu/pub/gsea/annotations/RN_U34.chip" # a 73 KB chip file
		con2 <- url(url.that.should.exist)
		tmp <- NULL
		try(tmp <- readLines(con2, 1), silent=TRUE)
		close(con2)
		if( is.null(tmp) ) {
			success <- FALSE
			cat("Connection to Broad FTP cannot be made.\n")
		}
		else {
			success <- FALSE
			cat("Your URL does not exist. (a connection to Broad can be made).\n")
		}
	}
	return( success )
}
