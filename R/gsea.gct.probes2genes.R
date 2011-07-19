#
# Take a table of expression data (at the probe-level), and convert it to a table of expression data
# at the unique gene-symbol-level.
# This chooses the best probe for each gene by looking at the stats vector which should be equivalent
# to the contents of a rnk file (todo: allow a rnk file to be specified)
# There's 2 ways of chooseing the best probe - Mark's way (Best_probe), or the default GSEA way
# which is called "Max_probe".
#
# Parameters:
#	data: a data.frame of expression data
#	stats: a vector of statistics
#	chip: either the filename of a chip file, or the result of import.gsea.chip
#	mode: Best_probe or Max_probe
#
# Value:
#	a data.frame of expression data with 1 row per unique gene symbol, and gene symbols as the rownames
#	see also gsea.genesets2gct.
#
# Mark Cowley, 2009-07-27
#


##' Take a table of expression data (at the probe-level), and convert it to a
##' table of expression data at the unique gene-symbol-level. This chooses the
##' best probe for each gene by looking at the stats vector which should be
##' equivalent to the contents of a rnk file (todo: allow a rnk file to be
##' specified) There's 2 ways of chooseing the best probe - Mark's way
##' (Best_probe), or the default GSEA way which is called "Max_probe".
##' 
##' @param data a data.frame of expression data
##' @param stats a vector of statistics
##' @param chip either the filename of a chip file, or the result of
##'   import.gsea.chip
##' @param mode Best_probe or Max_probe
##' @return a data.frame of expression data with 1 row per unique gene symbol,
##'   and gene symbols as the rownames see also gsea.genesets2gct.
##' @author Mark Cowley, 2009-07-27
##' @export
gsea.gct.probes2genes <- function(data, stats, chip, mode=c("Best_probe", "Max_probe", "Median_of_probes")) {
	mode <- mode[1]
	
	if( is.null(names(stats)) ) {
		stop("you must name your stats vector")
	}
	if( is.null(rownames(data)) ) {
		stop("you must name have rownames in your data")
	}
	
	if( !is.matrix.like(chip) && is.file(chip) ) {
		chip <- import.gsea.chip(chip)
	}

	# Cleanup the chip object to only consider the subset of probes that point to valid genes
	# remove the '---' genes
	chip <- chip[!is.na(chip$Gene.Symbol), ]
	chip <- chip[!chip$Gene.Symbol %in% c("---", ""), ]
	# sometimes the Gene.Symbol is REALLY long...
	chip <- chip[nchar(chip$Gene.Symbol) < 256, ]
	
	#
	# keep only the probes that are mentioned in names(stats) and rownames(data) and chip[,1]
	# -- this removes rows that have no stats which can be very useful.
	#
	probes <- intersectN(rownames(data), names(stats), chip[,1])
	tmp <- data.frame(
		data[match(probes, rownames(data)),], 
		stats=stats[match(probes, names(stats))], 
		chip[match(probes, chip[,1]), ], stringsAsFactors=FALSE)

	if( mode == "Best_probe" ) {
		o <- order(abs(tmp$stats), decreasing=TRUE)
	}
	else if( mode == "Max_probe" ) {
		o <- order(tmp$stats, decreasing=TRUE)
	}
	else if( mode == "Median_of_probes" ) {
		o <- order(tmp$stats, decreasing=TRUE)
	}
	tmp <- tmp[o, ]

	
	# collapse to one row per gene.
	rows <- as.list(vector2hashtable(tmp$Gene.Symbol))
	if( mode == "Best_probe" || mode == "Max_probe" ) {
		rows <- sapply(rows, "[", 1)
	}
	else if( mode == "Median_of_probes" ) {
		.get.median <- function(x) {
			idx <- c(1:length(x))
			idx <- floor(median(idx))
			x[idx]
		}
		rows <- sapply(rows, "[", .get.median)
	}
	tmp <- tmp[rows, ]
	
	# only keep the columns that refer to expression data.
	res <- tmp[, colnames(data)]
	rownames(res) <- tmp$Gene.Symbol
	
	return( res )
}


gsea.collapse.probes <- function(gct.file, rnk.file, chip.file, gct.file.out=NULL, rnk.file.out=NULL, mode=c("Best_probe", "Max_probe", "Median_of_probes")) {
	warning("code is untested. 2010-04-19.\n")
	mode <- mode[1]

	gct <- import.gsea.gct(gct.file)
	if( !missing(rnk.file) ) {
		rnk <- import.gsea.rnk(rnk.file) # named vector of statistics
	}
	else {
		rnk <- apply(gct[,3:ncol(gct)], 1, mean, na.rm=TRUE)
		names(rnk) <- gct[,1]
	}
	chip <- import.gsea.chip(chip.file)
	
	# Cleanup the chip object to only consider the subset of probes that point to valid genes
	# remove the '---' genes
	chip <- chip[!is.na(chip$Gene.Symbol), ]
	chip <- chip[!chip$Gene.Symbol %in% c("---", ""), ]
	# sometimes the Gene.Symbol is REALLY long...
	chip <- chip[nchar(chip$Gene.Symbol) < 256, ]
	
	#
	# keep only the probes that are mentioned in names(rnk) and rownames(gct) and chip[,1]
	# -- this removes rows that have no rnk which can be very useful.
	#
	probes <- intersect(rownames(gct), intersect(names(rnk), chip[,1]))
	tmp <- data.frame(
		gct[match(probes, rownames(gct)),], 
		rnk=rnk[match(probes, names(rnk))], 
		chip[match(probes, chip[,1]), ]
	)

	if( mode == "Best_probe" ) {
		o <- order(abs(tmp$rnk), decreasing=TRUE)
	}
	else if( mode == "Max_probe" ) {
		o <- order(tmp$rnk, decreasing=TRUE)
	}
	else if( mode == "Median_of_probes" ) {
		o <- order(tmp$rnk, decreasing=TRUE)
	}
	tmp <- tmp[o, ]

	#
	# choose one row per gene.
	#
	rows <- as.list(vector2hashtable(tmp$Gene.Symbol))
	if( mode == "Best_probe" || mode == "Max_probe" ) {
		rows <- sapply(rows, "[", 1)
	}
	else if( mode == "Median_of_probes" ) {
		.get.median <- function(x) {
			idx <- c(1:length(x))
			idx <- floor(median(idx))
			x[idx]
		}
		rows <- sapply(rows, "[", .get.median)
	}
	
	#
	# re-create the gct & rnk files
	#
	gct.out <- tmp[rows, colnames(gct)]
	rnk.out <- tmp$rnk; names(rnk.out) <- rownames(tmp)
	
	#
	# map identifiers to gene symbols, or leave as probeset id's?
	#
	
	if( is.null(gct.file.out) ) gct.file.out <- gct.file
	if( is.null(rnk.file.out) ) rnk.file.out <- rnk.file
	
}
