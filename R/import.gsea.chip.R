#' Import a GSEA chip file.
#' This can look inside a gsea output dir, find the path to the chip file that
#' was used, and try to open it.
#' 
#' Missing values in chip files.\cr
#' There seems to be little convention re missing values in chip files. It's important
#' to map all potential missing values to \code{NA}'s, since lots of code in this package 
#' assumes that the symbol is either valid, or \code{NA}. Set \code{na.strings} accordingly.
#' 
#' Gene Symbol capitalisation.\cr
#' GSEA likes to capitalise all gene symbols, so in GSEA result files, the symbols will be
#' all caps. This can cause logical problems when using the symbol as a key. 
#' So, even though some human gene symbols are not supposed to be in caps, eg C1orf34,
#' we recommend that you set \code{capitalise.symbol=TRUE}.
#' 
#' @param x the path or url to a chip file, or the path to a gsea directory
#' @param na.strings a character vector of values which should be mapped to NA. see details
#' @param capitalise.symbol logical: capitalise all gene symbols upon import? see details
#' @param verbose logical: if TRUE, then write out verbose messages, including download.file
#'    verbose output.
#' @author Mark Cowley, 2009-08-07
#' @export
import.gsea.chip <- function(x, na.strings=c("---", "NA", ""), capitalise.symbol=TRUE, verbose=TRUE ) {
	if( is.gsea.dir(x) ) {
		rpt <- import.gsea.rpt(x)
		import.gsea.chip( rpt$chip, na.strings=na.strings, capitalise.symbol=capitalise.symbol )
	}
	else if( is.url(x) ) {
		f <- tempfile()
		if( verbose ) cat("Downloading .chip file.\n")
		success <- download.file(x, f, quiet=!verbose)
		success == 0 || stop(sprintf("Could not download chip file at: '%d'", x))
		import.gsea.chip(f, na.strings=na.strings, capitalise.symbol=capitalise.symbol)
	}
	else if( file.exists(x) ) {
		chip <- read.delim(x, stringsAsFactors=FALSE)
		for(na in na.strings)
			chip[chip == na] <- NA
		chip <- chip[, 1:match("Gene.Title", colnames(chip))]
		if( capitalise.symbol ) {
			chip$Gene.Symbol <- toupper(chip$Gene.Symbol)
		}
			
		if( is.numeric(chip$Probe.Set.ID) ) {
			# if( verbose ) message("Probe ID's are numeric - converting them to text.\n")
			chip$Probe.Set.ID <- as.character(chip$Probe.Set.ID)
		}

		if( any(nchar(chip$Gene.Symbol)>255) ) {
			idx <- which(nchar(chip$Gene.Symbol)>255)
			tmp <- substring(chip$Gene.Symbol[idx], 1, 255)
			# warning: this will frequently truncate the last gene symbol in the string.
			# todo: cleanup tmp by nibbling from the end to ensure the strings end on a proper gene symbol
			chip$Gene.Symbol[idx] <- tmp
		}
		chip
	} 
	else {
		stop("Can't locate the chip file:", x)
	}
}
# CHANGELOG
# 2009-08-07: v1
# 2011-02-22: if chip[,1] is.numeric, then make it character (this is the case for Affy Gene & Exon ST)
