#' Import a cls file
#'
#' Import a GenePattern cls file into R. CLS files can either be categorical, or continuous & 
#' both can be handled by this code. 
#' 
#' The cls file spec says for catgegorical cls files, values should start from 0. 
#' However some tools (eg SubMap) require them to start from 1. \code{enforce.zero.start} 
#' allows control over these situations.
#'
#' @param file  the path to a valid cls file.
#' @param as.factor logical: return a character, or factor?
#' @param enforce.zero.start logical: check that the lowest class is 0. See details.
#' @return if cls file is continuous, either a character vector, or factor, one element per sample,
#' of if continuous, a numeric is returned.
#' @author Mark Cowley, 2009-11-26
#' @seealso \code{\link{import.gsea.gct}} \code{\link{import.gsea.clm}} \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#cls}
#' 
import.gsea.cls <- function(file, as.factor=TRUE, enforce.zero.start=TRUE) {
	raw <- readLines(file, warn=FALSE)
	raw <- trim(raw)

	if( length(raw) <= 2 ) {
		stop("cls file has too few rows.\n")
	}

	if( raw[1] == "#numeric" ) {
		raw <- raw[2:length(raw)]
		raw <- raw[nchar(raw) > 0] # skip blank rows.
		N <- length(raw)
		stopifnot( even(N) )
		
		names <- raw[seq(1,N,2)]
		names <- sub("^#[ \t]*", "", names)
		
		values <- raw[seq(2,N,2)]
		values <- strsplit(values, "[ \t]")
		values <- lapply(values, as.numeric)
		values <- t( as.data.frame(values) )
		rownames(values) <- names
		colnames(values) <- NULL

		return( values )
	}
	else { # categorical cls file, which contains only one label.
		# line 1
		tmp <- strsplit(raw[1], "[ \t]")[[1]]
		Nsamples <- tmp[1]
		Nclasses <- tmp[2]
	
		# line 2
		stopifnot(substr(raw[2], 1, 1) == "#")
		raw[2] <- sub("#[ \t]*", "", raw[2])
		classes <- strsplit(raw[2], "[ \t]")[[1]]
		stopifnot(length(classes) == Nclasses)
	
		# line 3. 
		tmp <- strsplit(raw[3], "[ \t]")[[1]]
		values <- as.numeric(tmp)
		if( enforce.zero.start && min(values) != 0 )
			stop("Minimum value in cls file must be zero. Set enforce.zero.start=FALSE to override.\n")
		if( min(values) == 0 ) values <- values + 1
		stopifnot(length(values) == Nsamples)
		stopifnot( max(values) == Nclasses ) # NB i've added 1 already if values started at 0
		res <- classes[values]
		if( as.factor )
			res <- factor(res, levels=classes)
	
		# tmp <- strsplit(raw[3], "\t")[[1]]
		# values <- tmp
		# stopifnot(length(values) == Nsamples)
		# unique.values <- unique(values)
		# stopifnot(length(unique.values) == Nclasses)
		# res <- factor(values, classes=unique.values)
		# classes(res) <- classes
		# 
		# if( !as.factor )
		# 	res <- as.character(res)

		res
	}
}
# CHANGELOG
# 2010-07-06: update that allows minimum value to be 1 (for using Submap)
#