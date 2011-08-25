#' Is the argument a GSEA zip file?
#'
#' @param x a vector of file paths
#' @return a logical vector, same length as \code{x}
#' @author Mark Cowley, 2011-08-24
#' @export
#' @examples
#' f <- "my.zip"
#' is.gsea.zip(f)
#' is.gsea.zip(c(f,f,f,"/home"))
is.gsea.zip <- function(x) {
	is.file(x) & grepl("zip$", x, ignore.case=TRUE)
}
