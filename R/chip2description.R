#' Convert a GSEA chip file or object into a vector of descriptions.
#' @param chip.file the path to a chip file; or:
#' @param chip an imported chip object (see \code{\link{import.gsea.chip}})
#' @param sep the field separator between the symbol and description
#' @param probes an optional vector of probe names to limit the result to. Default is \code{NULL}
#'   whereby all probes will be converted
#' @param genoloc optional vector of genomic locations, \code{NA} and \code{NULL} allowed, MUST be in
#'    identical order to the rows of the chip
#' @return a character vector of description. If genoloc is provided, then eg \dQuote{|@@chr6:1950428-1950681|}
#'  will be appended.
#' @author Mark Cowley
#' @export
chip2description <- function(chip.file=NULL, chip=NULL, sep=": ", probes=NULL, genoloc=NULL) {
	if( !is.null(chip.file) && length(chip.file) == 1 && (file.exists(chip.file) || is.url(chip.file)) ) {
		chip <- import.gsea.chip(chip.file)
	}
	else if( !is.null(chip.file) && is.data.frame(chip.file) ) {
		stop("You passed chip object where a chip filename was expected.\n")
	}
	description <- apply(chip[,2:3], 1, paste, collapse=sep)
	names(description) <- as.character(chip[,1])
	description[is.na(chip[,2]) & is.na(chip[,3])] <- ""
	description[!is.na(chip[,2]) & is.na(chip[,3])] <- chip[!is.na(chip[,2]) & is.na(chip[,3]), 2]
	description[is.na(chip[,2]) & !is.na(chip[,3])] <- chip[is.na(chip[,2]) & !is.na(chip[,3]), 3]
	
	if( !is.null(probes) ) {
		probes <- as.character(probes)
		# limit & reorder the result to just these probes, also adding blank descriptions
		# for probes that are not in the chip file.
		if( any(! probes %in% chip[,1]) ) {
			ids <- probes[which(!probes %in% chip[,1])]
			cat(sprintf("WARNING: %d probes in the data were not referenced in the chip file.\n", length(ids)))
			tmp <- rep("", length(ids))
			names(tmp) <- ids
			description <- c(description, tmp)
			# description <- description[probes]
		}
		if( !is.null(genoloc) ) genoloc <- genoloc[match(probes, names(description))]
		description <- description[match(probes, names(description))]
	}

	#
	# if genomic locations are provided as a character vector in ucsc format, 
	#   eg: character("chr6:1950428-1950681", ...)
	# then append " |@chr6:1950428-1950681|" to the end of the genes' description.
	# >= 0 genloc's are OK, but >1 should be comma separated, eg "chr6:1950428-1950681,chr6:2304548-2304574"
	#
	if( !is.null(genoloc) ) {
		length(genoloc) == nrow(chip) || stop("length genoloc must match number of rows in chip, after optional probe filtering")
		genoloc[is.na(genoloc) | is.null(genoloc)] <- ""
		genoloc <- trim(genoloc)
		genoloc <- sub("chrom", "chr", genoloc)	# probably never used
		genoloc <- gsub(", +", ",", genoloc)	# multiple genolocs per probe are OK, but must be csv with no spaces allowed.
		genoloc <- paste("|@",genoloc,"|", sep="")
		genoloc[genoloc == "|@|"] <- ""
		description <- paste(description, genoloc)
		description <- trim(description) # trim the " " which will be added to end if no genoloc provided.
	}

	return( description )
}
