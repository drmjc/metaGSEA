#' Function to export a GenePattern gct file.
#' 
#' \code{data} can be either:\cr
#' a \code{matrix} or \code{data.frame} of numeric values, 
#' in which case you can provide an optional \code{description} or \code{chip} 
#' object to populate the \dQuote{Description} column;
#' or a \code{data.frame} containing \sQuote{Name} 
#' and \dQuote{Description} columns (see \code{\link{import.gsea.gct}})
#' 
#' @param data a \code{matrix} or \code{data.frame} WITH rownames of all numeric
#'   values, or a \code{data.frame} from running \code{\link{import.gsea.gct}}; 
#'  ie first 2 columns are \dQuote{Name}, \dQuote{Description}.
#' @param description a vector of annotations. IF it has names, then these
#'   names must match the rownames of data, in which case we make sure they're
#'   in the same order. If no names, then we assume that they're in the same
#'   order. If \code{NULL}, then you must either provide a \code{chip}, or have a
#'   \dQuote{Description} column in the data.
#' @param file the output file name
#' @param chip instead of specifying description, you can specify a GSEA chip
#'   object, and a description will be made for you. This overrides
#'   \code{description}. Default is \code{NULL} to ignore.
#' @param round the number of digits to round the numbers to - default=4
#' @param version The GCT file version, to go in the first line.
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
			if( anyDuplicated(data$Name) != 0 )
				stop("Values in 'data$Name' MUST be unique as per the gct file format definition.\n")
			rownames(data) <- data$Name
			description <- data$Description; names(description) <- rownames(data)
			data <- data[,3:ncol(data), drop=FALSE]
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
# CHANGELOG
# 2011-11-09
# - Allow single sample GCT files to work
# - named arguments across my workspace and GenePattern modules
# - @TODO: move the file argument to #2


# 
# export.broad.gct <- function(data, description, file, round=3, version="#1.2", missing="", ...) {
# 	export.gsea.gct(data=data, description=description, file=file, round=round, version=version, missing=missing, ...)
# }
