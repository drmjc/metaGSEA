#' import.genepattern.odf
#' 
#' Import a GenePattern ODF file. This has been tested using ODF created by
#' LimmaGP, which in turn was modelled on that from ComparativeMarkerSelection.
#' Note that all the header info is stripped out & lost upon import,
#' except for the column names.
#'
#' @param f the path to a single ODF file
#' @return a \code{data.frame} representation of the ODF file.
#' @author Mark Cowley, 2011-11-22
#' @export
#' @importFrom stringr str_replace
#' @examples
#' \dontrun{
#' f <- "yellow_exp_LeanConvsObCon.odf"
#' a <- import.genepattern.odf(f)
#' head(a)
#' }
import.genepattern.odf <- function(f) {
	header <- readLines(f, 100)
	grepl("ODF 1.0", header[1]) || stop("Doesn't look like an ODF 1.0 file")
	any(grepl("HeaderLines", header)) || stop("HeaderLines not found in ODF header")
	any(grepl("DataLines", header)) || stop("DataLines not found in ODF header")
	any(grepl("COLUMN_NAMES", header)) || stop("COLUMN_NAMES not found in ODF header")
	
	HeaderLines <- grep("HeaderLines", header, value=TRUE)
	HeaderLines <- as.numeric(str_replace(HeaderLines, ".*=", ""))
	!is.na(HeaderLines) || stop("HeaderLines=?? not a number")
	skip <- HeaderLines + 1 # skip the ODF 1.0 row + the num rows specified my HeaderLines (does not include the DataLines row)
	
	col.names <- grep("COLUMN_NAMES", header, value=TRUE)
	col.names <- str_replace(col.names, "COLUMN_NAMES:", "")
	col.names <- strsplit(col.names, "\t")[[1]]
	
	grep("DataLines", header) == (skip+1) || stop("DataLines not in the expected spot")
	
	DataLines <- grep("DataLines", header, value=TRUE)
	DataLines <- as.numeric(str_replace(DataLines, ".*=", ""))
	!is.na(DataLines) || stop("DataLines=?? not a number")
	skip <- skip + 1 # skip the DataLines row
	
	res <- read.delim(f, header=FALSE, skip=skip)
	nrow(res) == DataLines || stop(sprintf("Found %d lines of data, expected %d", nrow(res), DataLines))
	colnames(res) <- col.names
	
	return( res )
}
# CHANGELOG
# 2011-11-29:
# - from limmaGP, the HeaderLines doesn't include the DataLines row. CMSV likes LimmaGP output, but still, need to check the CMS outputted ODF's (@TODO)