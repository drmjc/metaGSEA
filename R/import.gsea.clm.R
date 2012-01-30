#' Import a clm file
#'
#' Import a GenePattern clm file. CLM files are often used to annotate GCT files.
#'
#' @param f character: the path to a valid clm file.
#' @author Mark Cowley, 2009-12-07
#' @seealso \code{\link{import.gsea.gct}} \code{\link{import.gsea.cls}} 
#' \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#clm}
#' @return a \code{data.frame} with 3 columns: \sQuote{CEL}, \sQuote{Sample}, \sQuote{Class}
#' @export
import.gsea.clm <- function(f) {
	suppressWarnings(
		clm <- read.delim(f, header=FALSE, stringsAsFactors=FALSE, colClasses = "character", comment.char="")
	) # suppressWarnings implemented to avoid warnings if the last line is missing a newline. This could suppress genuine warnings though.
	if(ncol(clm) != 3)
		stop("clm files should have just three columns. watch out for extra tabs.\n")
	if( any(is.na(clm)) )
		stop("There should be no missing fields in a clm file.\n")
	if( length(grep(" +$", clm[,1])) > 0 )
		warning("Some of the sample names in column 1 of the CLM file have a trailing space. It's safest if you trim these trailing spaces.\n")

	colnames(clm) <- c("CEL", "Sample", "Class")
	if( length(unique(clm$Sample)) != length(clm$Sample) ) {
		clm$Sample <- make.unique(clm$Sample)
		warning("Your sample names (column 2) were not unique! I've made them unique.\n")
	}
	if( !all(clm$Class == make.names(clm$Class, unique=FALSE)) ) {
		warning("The class names in the third column in clm file must start with a letter, and can only contain letters, numbers, '.' and '_'\nAll suspect characters have been changed to '.'\n")
		clm$Class <- make.names(clm$Class, unique=FALSE)
	}

	return( clm )
}
# 2012-01-24
# - added check for trailing spaces in clm[,1]
