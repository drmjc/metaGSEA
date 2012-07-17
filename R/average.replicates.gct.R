#' Average replicate columns from within a gct file
#' Average replicates in a gct file, based on the
#' classes in \code{classes}
#'
#' @section TODO:
#' when I have a GCT file, make this use the S4 method
#' 
#' @param gct.file the path to a gct file.
#' @param classes a character or factor: the sample classes, 1 per column in x, with as many
#' levels as there are unique samples
#' @return A data.frame from importing the gct file, and then restricting the data.
#' @author Mark Cowley, 2011-09-01
#' @export
#' @importMethodsFrom microarrays average.replicates
average.replicates.gct <- function(gct.file, classes) {

	x <- import.gsea.gct(gct.file)
	lhs <- x[,1:2]
	x <- x[,3:ncol(x)]
	
	res <- average.replicates(x, classes)
	res <- data.frame(lhs, res)

	res
}
