# Convert a clm file into a cls file.
#
# Parameters:
#	clm.file: the path to a clm file
#	cls.file: if NULL, the cls file will be in the same dir as the clm, and will have the suffix cls.
#
# Mark Cowley, 2009-12-21


##' Convert a clm file into a cls file.
##' 
##' @param clm.file the path to a clm file
##' @param cls.file if NULL, the cls file will be in the same dir as the clm,
##'   and will have the suffix cls.
##' @author Mark Cowley, 2009-12-21
##' @export
gsea.convert.clm2cls <- function(clm.file, cls.file=NULL) {
	if( is.null(cls.file)) 
		cls.file <- sub("clm", "cls", clm.file)
		
	tmp <- import.gsea.clm(clm.file)
	export.broad.cls(tmp[,3], cls.file, continuous=FALSE)
	
	invisible(cls.file)
}
