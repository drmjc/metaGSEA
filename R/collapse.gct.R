#' Collapse a GCT file
#' 
#' It's common for microarrays to have multiple probes per gene. They tend to
#' represent different isoforms.
#' Most geneset testing is done at the gene symbol level & ignores isoforms,
#' so you need to choose 1 probe
#' for each gene. How? 2 common approaches are to take the most abundant
#' probe, or the most variable probe, considered
#' across the cohort.
#' I quite like doing t-stats on each gene & selecting the best performing
#' probe - ie the one with the largest t-stat in
#' either direction. Why? On the Affy 133+2 array, there can be lots of poor
#' probes for each gene. If 5 probes for a gene
#' have these t-stats: 1.2, 0.9, 0.1, -0.1, -10; then IMO, the one that scored
#' -10 is the best probe, since it had a really
#' strong t-stat score. thus \code{method="maxabs"} combined with a rnk.file
#' 
#' @param gct.file the path to a gct file
#' @param chip.file the path to a chip file
#' @param gct.outfile the path to the gct output file
#' @param rnk.file [optional] path to a rnk file (eg a t-statistic for each
#'   probe, where you want to select best probe from this score) NB currently
#'   UNUSED 
#' @param method \dQuote{mean}, \dQuote{median} select the probe with highest 
#'   average/median level, or \dQuote{var}: select the probe with highest variance 
#'   across samples; \dQuote{maxabs} select the probe with the large absolute score in the rnk
#'   file (see details).
#' @param reverse [default=FALSE] reverse the ordering selected by method arg.
#'   so instead of most variable, it would be least variable.
#' @param filter Filter out (ie exclude) those probes that don't have a gene
#'   symbol (as determined by probes that have a gene symbol of \code{NA}, \dQuote{},
#'   \dQuote{---}, or \dQuote{NA}.)
#' @return A gct file is created with 1 row per gene symbol & now the \sQuote{probe ids}
#'  in column 1 are actually gene symbols.
#' @author Mark Cowley, 2011-02-27
#' @export
collapse.gct.file <- function(gct.file, chip.file, gct.outfile, rnk.file=NULL, method=c("var", "mean", "median"), reverse=FALSE, filter=FALSE) {
	gct <- import.gsea.gct(gct.file)
	chip <- import.gsea.chip(chip.file)
	rnk <- NULL
	if( !is.null(rnk.file) ) rnk <- import.gsea.rnk(rnk.file)
	res <- collapse.gct(gct, chip, rnk, method=method, reverse=reverse, filter=filter)
	export.gsea.gct(res, file=gct.outfile)
}

#' Collapse a GCT object
#' 
#' It's common for microarrays to have multiple probes per gene. They tend to
#' represent different isoforms.
#' Most geneset testing is done at the gene symbol level & ignores isoforms,
#' so you need to choose 1 probe
#' for each gene.
#' How?
#' 2 common approaches are to take the most abundant probe, or the most
#' variable probe, considered
#' across the cohort.
#' I quite like doing t-stats on each gene & selecting the best performing
#' probe - ie the one with the largest t-stat in
#' either direction. Why? On the Affy 133+2 array, there can be lots of poor
#' probes for each gene. If 5 probes for a gene
#' have these t-stats: 1.2, 0.9, 0.1, -0.1, -10; then IMO, the one that scored
#' -10 is the best probe, since it had a really
#' strong t-stat score. thus method="maxabs" combined with a rnk.file (NB NOT
#' implemented right now.)
#' 
#' @param gct a GCT object
#' @param chip a CHIP object
#' @param rnk [optional] path to a rnk file (eg a t-statistic for each probe,
#'   where you want to select best probe from this score) NB currently UNUSED
#' @param method \dQuote{mean}, \dQuote{median} select the probe with highest 
#'   average/median level, or \dQuote{var}: select the probe with highest variance 
#'   across samples; \dQuote{maxabs} select the probe with the large absolute score in the rnk
#'   file (see details).
#' @param reverse [default=FALSE] reverse the ordering selected by method arg.
#'   so instead of most variable, it would be least variable.
#' @param filter Filter out (ie exclude) those probes that don't have a gene
#'   symbol (as determined by probes that have a gene symbol of \code{NA}, \dQuote{},
#'   \dQuote{---}, or \dQuote{NA}.)
#' @return A gct object with 1 row per gene symbol & now the \sQuote{probe ids} in
#'   column 1 are actually gene symbols.
#' @author Mark Cowley, 2011-02-27
#' @export
collapse.gct <- function(gct, chip, rnk=NULL, method=c("var", "mean", "median"), reverse=FALSE, filter=FALSE) {
	if( !all(gct[,1] %in% chip[,1]) ) {
		stop("Some ID's in gct are not in the chip file.")
	}
	else {
		chip <- chip[match(gct[,1], chip[,1]), ]
	}
		
	#
	# reorder the rows of the gct so that better probes are higher up that worse probes.
	# (better is set according to the given 'method')
	#
	gct <- reorder.gct(gct, method=method, reverse=reverse)
	if( reverse ) gct <- gct[nrow(gct):1, ]
	chip <- chip[match(gct[,1], chip[,1]), ]
	gct$ORDER <- 1:nrow(gct) # used to integrate the rows that do and do NOT have probe symbols.
	
	# split the data into those probes that do & do not have Gene Symbols
	probes.no.sym <- chip[,2] == "NA" | chip[,2] == "---" | chip[,2] == "" | is.na(chip[,2])
	gct.hasSymbols <- gct[!probes.no.sym, ]
	gct.noSymbols  <- gct[probes.no.sym, ]
	chip.hasSymbols <- chip[!probes.no.sym, ]
	chip.noSymbols  <- chip[probes.no.sym, ]

	# for the probes that DO have a valid gene symbol, choose the best probe.
	uSymbols <- unique(chip.hasSymbols[,2])
	m <- match(uSymbols, chip.hasSymbols[,2])
	res <- data.frame( chip.hasSymbols[m, 2:3], 
					   gct.hasSymbols[m, 3:ncol(gct.hasSymbols)], 
					   check.names=FALSE )
	colnames(res)[1:2] <- colnames(gct)[1:2]
	rownames(res) <- res[,1]
	
	# for the probes that DO NOT have a valid gene symbol, use the same information.
	if( !filter && nrow(gct.noSymbols) > 0 ) {
		# Reannotate the first 2 columns using the chip file
		gct.noSymbols[,2] <- chip2description(chip=chip.noSymbols)
		res <- rbind(res, gct.noSymbols)
		res <- res[order(res$ORDER), ]
	}

	res$ORDER <- NULL
	
	res
}
