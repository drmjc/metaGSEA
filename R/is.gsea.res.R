#' is.gsea.gct.file
#'
#' @param f the path to a file
#' @return \code{TRUE} if \code{f} is a GCT file, \code{FALSE} otherwise
#' @author Mark Cowley, 2011-10-21
#' @export
is.gsea.gct.file <- function(f) {
	grepl("gct$", f, ignore.case=TRUE)
}

#' is.gsea.res.file
#'
#' @param f the path to a file
#' @return \code{TRUE} if \code{f} is a RES file, \code{FALSE} otherwise
#' @author Mark Cowley, 2011-10-21
#' @export
is.gsea.res.file <- function(f) {
	grepl("res$", f, ignore.case=TRUE)
}
