#' is the argument a valid GSEA directory?
#'
#' @param x path to a directory
#' @return \code{TRUE} if \code{x} is the top-level GSEA result directory, \code{FALSE} otherwise
#' @author Mark Cowley, 2011-08-04
#' @export
is.gsea.dir <- function(x) {
	!is.null(x) & !is.na(x) & file.exists(x) & is.dir(x) & file.exists(file.path(x, "index.html"))
}
# is.gsea.dir.dir <- function(x) {
# 	if( file.exists(x) && is.dir(x) && !is.gsea.dir(x) ) {
# 		dirs <- dir(x, full.path=TRUE)
# 		dirs <- dirs[is.dir(dirs)]
# 		any(sapply(dirs, is.gsea.dir))
# 	}
# 	FALSE
# }
# 