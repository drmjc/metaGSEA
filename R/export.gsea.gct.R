#' Function to export a GenePattern gct file.
#' 
#' \code{data} can be either a \code{matrix} or \code{data.frame} of numbers, 
#' in which case description and or chip should be provided to annotate the rows; 
#' or a \code{data.frame} containing a \sQuote{Name} 
#' and \dQuote{Description} columns (see \code{\link{import.gsea.gct}})
#' 
#' @param data a data.frame WITH rownames, and optionally can be from
#'   import.gsea.gct; ie first 2 columns are Name, Description
#' @param description a vector of annotations. IF it has names, then these
#'   names must match the rownames of data, in which case we make sure they're
#'   in the same order. If no names, then we assume that they're in the same
#'   order. If NULL, then you must either provide a chip, or have a
#'   Description column in the data.
#' @param file the output file name
#' @param chip instead of specifying description, you can specify a GSEA chip
#'   object, and a description will be made for you. This overrides
#'   description. Default is NULL to ignore.
#' @param round the number of digits to round the numbers to - defaults to 3
#' @param version The first line of the gct file will have this tag.
#' @param missing the string to use for missing data (in the expression data).
#' @param \dots	Currently unused.
#' @author Mark Cowley, 2008-08-07
#' @examples
#' \dontrun{
#' my.gct <- import.gsea.gct("./my.gct", "./my.cls")
#' export.gsea.gct(data=my.gct, file="my2.gct")
#' }
#' @export
export.gsea.gct <- function(data, description=NULL, file=NULL, chip=NULL, round=4, version="#1.2", missing="", ...) {
	if( is.null(description) && is.null(chip) ) {
		if(!"Description" %in% colnames(data)) {
			description <- rep("", nrow(data))
			# stop("You did not provide a description or chip, nor is there a column named Description in your data.\n")
		}
		else {
			# data looks like its a GSEA gct object
			if( length(data$Name) != length(unique(data$Name)) )
				stop("Values in data$Name MUST be unique as per the gct file format definition.\n")
			rownames(data) <- data$Name
			description <- data$Description; names(description) <- rownames(data)
			data <- data[,3:ncol(data)]
		}
		export.gsea.gct(data, description=description, file=file, chip=NULL, round=round, version=version, missing=missing, ...)
	}
	else if( !is.null(chip) ) {
		description <- chip2description(chip=chip, probes=rownames(data))
		export.gsea.gct(data, description=description, file=file, chip=NULL, round=round, version=version, missing=missing, ...)
	}
	else {
		data <- round(data, round)
	
		na.count <- rowSums(is.na(data))
		if( any(na.count == ncol(data)) ) {
			warning("Some rows were all NA - these will be excluded.\n")
			data <- data[,na.count < ncol(data)]
			description <- description[na.count < ncol(data)]
		}
		
		if( is.null(names(description)) && length(description) != nrow(data) ) {
			stop("description has no names, and is a different length to nrow(data). It's best if you name the description elements using the probeset id.\n")
		}
		else if( ! is.null(names(description)) ) {
			if( all(rownames(data) %in% names(description)) ) {
				description <- description[rownames(data)]
			}
			else {
				stop(sprintf("%d probes in data were not found in the description vector. Make sure that all probes in the data have an element in the description vector.", sum(!rownames(data) %in% names(description))))
			}
		}

		gct <- data.frame(Name=rownames(data), Description=description, data, stringsAsFactors=FALSE, check.names=FALSE)
	
		OUT <- file(file, "w")
		writeLines(version, OUT)
		writeLines(paste(nrow(data), ncol(data), sep="\t"), OUT)
		write.delim(gct, OUT, na=as.character(missing))
		close(OUT)
		
		invisible( gct )
	}
}
# 
# export.broad.gct <- function(data, description, file, round=3, version="#1.2", missing="", ...) {
# 	export.gsea.gct(data=data, description=description, file=file, round=round, version=version, missing=missing, ...)
# }
