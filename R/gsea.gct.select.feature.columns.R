#' Restrict the number of columns in a gct & cls file.
#' 
#' You can do this either by indexing the columns ([1:n]), or by specifying the
#' names of the
#' classes. The order of columns, or classes is respected.
#' @param gct.file the path to a GCT file
#' @param cls.file the path to a CLS file
#' @param gct.out the path to a GCT output file
#' @param cls.out the path to a CLS output file
#' @param columns the column indices to include, or \code{NULL}
#' @param classes the classes to include, or \code{NULL}
#' @param cls.start.value most cls files start from 0, but that can be changed if necessary
#' @return none. creates a gct and cls file containing the specified subset of classes or columns.
#' @author Mark Cowley, 2010-07-06
#' @export
gsea.gct.select.features.columns <- function(gct.file, cls.file, gct.out, cls.out, columns=NULL, classes=NULL, cls.start.value=0) {
	if( is.null(columns) && is.null(classes) )
		stop("Must specify either columns, or classes.\n")

	gct <- import.gsea.gct(gct.file)
	cls <- import.gsea.cls(cls.file, enforce.zero.start=FALSE)
	if( length(cls) != ncol(gct)-2 ) {
		stop("number of samples in gct and cls do not match.\n")
	}
	
	if( !is.null(classes) ) {
		if( any(! classes %in% levels(cls) ) ) {
			stop( sprintf(
					"These classes were not found in the cls file: %s", 
					tocsv( setdiff(classes, levels(cls)) )
			) )
		}
		# to respect the order of classes, need to do this in a loop.
		columns <- which(cls %in% classes[1])
		for(class in classes[2:length(classes)]) {
			columns <- c(columns, which(cls %in% class))
		}
	}
	else {
		if( min(columns <= 0) || max(columns) > length(cls) )
			stop("columns must be in [1,n] where n is the number of samples in the cls file.\n")
	}
	
	gct <- gct[,c(1,2,columns+2)]
	cls <- factor(as.character(cls)[columns], levels=classes)
	export.gsea.gct(gct, file=gct.out)
	export.broad.cls(cls, file=cls.out, start.value=cls.start.value)
}
