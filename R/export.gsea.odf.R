#' Export a table of data to a GenePattern odf file
#' 
#' ODF files are tab delimited documents with fairly extensive header lines. Strongly
#' recommend you look at the documentation for odf files (see below), and to see it
#' in action, see \code{\link{export.gsea.odf.lmFit}}.
#' 
#' @note header lines.\cr
#' the header lines should not contain "ODF 1.0", "HeaderLines=X", though the code 
#' ignore them.
#' The only required field that you should consider providing is the
#' \dQuote{Model=DataSet}.
#' The \dQuote{DataLines=XYZ} header line is determined automatically.
#' 
#' @param data a \code{\link{data.frame}} or \code{\link{matrix}}
#' @param file the path to the output odf file
#' @param header.lines a character vector of header lines, or \code{NULL}
#' @param colclasses a character vector of column class names. if \code{NULL}, then
#'    they are aut-determined from the \code{colclass} of each column, and mapped
#'    to GenePattern column class names
#' @return none. writes a file
#' @author Mark Cowley, 2009-11-30
#' @seealso \code{\link{export.gsea.odf.lmFit}} \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#odf}
#' @export
#' 
export.gsea.odf <- function(data, file, header.lines=NULL, colclasses=NULL) {
	
	######################################
	# Header Lines checking
	# -- strip these values out, as they will be written later on.
	# this will leave any non-standard header lines
	#
	if( any(grepl("ODF 1.0", header.lines)) ) {
		header.lines <- header.lines[-grep("ODF 1.0", header.lines)]
	}
	if( any(grepl("HeaderLines=", header.lines)) ) {
		header.lines <- header.lines[-grep("HeaderLines=", header.lines)]
	}
	if( !any(grepl("Model=", header.lines)) ) {
		header.lines <- c(header.lines, "Model=Dataset")
	}
	if( !any(grepl("DataLines=", header.lines)) ) {
		header.lines <- c(header.lines, paste("DataLines", nrow(data), sep="="))
	}
	
	if( is.null(colnames(data)) ) {
		colnames(data) <- paste("V", 1:ncol(data), sep="")
	}
	# end header lines checking
	######################################

	
	######################################
	# ascertain the data type of each column. factors and logicals will be converted.
	#
	if( is.null(colclasses) ) {
		types <- colclasses(data)
		if( any(types == "factor") ) {
			for(i in which(types=="factor")) {
				data[,i] <- as.character(data[,i])
			}
		}
		if( any(types == "logical") ) {
			for(i in which(types=="logical")) {
				data[,i] <- as.numeric(data[,i])
			}
		}
		types <- colclasses(data)

		types[types=="character"] <- "String"
		types[types=="numeric"] <- "double"
		types[types=="integer"] <- "int"
		#
		# R calls doubles and integers within a data.frame just "numeric"...
		# do an integer test on the numerics.
		#
		row1 <- trim(as.character(unlist(data[1,])))
		if( any(is.int(row1)) ) { # ignore columns that aren't int's in the very first row.
			idx <- which(types== "double" & is.int(data[1,]))
			if( length(idx) > 0 ) {
				for(i in idx) {
					if( all(is.int(data[,i]) | is.na(data[,i])) )
						types[i] <- "int"
				}
			}
		}
	}
	else {
		if( any(!colclasses %in% c("int", "String", "double")) ) {
			stop("colclasses can only contain int, String or double's.\n")
		}
		types <- colclasses
	}
	# done checking for types.
	######################################

	### construct the header ##
	nHeaderLines <- 2 + length(header.lines)
	header <- c(
		"ODF 1.0",
		paste("HeaderLines=", nHeaderLines, sep=""),
		paste("COLUMN_NAMES:", paste(colnames(data), collapse="\t"), sep=""),
		paste("COLUMN_TYPES:", paste(types, collapse="\t"), sep=""),
		header.lines
		)
		
	### write out the file ###
	OUT <- file(file, "w")
	writeLines(header, OUT)
	write.delim(data, OUT, col.names=FALSE, row.names=FALSE)
	close(OUT)
}
