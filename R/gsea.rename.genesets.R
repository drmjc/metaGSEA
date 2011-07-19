# Rename the genesets from within a GSEA object.
#
# Rename genesets using a prefix and suffix, such that COWLEY_UP becomes: PREFIX_COWLEY_UP_SUFFIX
# Optionally, the rank, FDR and direction can also be included so that COWLEY_UP becomes:
#	1: PREFIX_COWLEY_UP_SUFFIX UP (0.002)
#
# Parameters:
#	x: a GSEA list, optionally with an edb. see import.gsea(..., edb=TRUE)
#	prefix: an optional prefix (such as the treatment, eg WTvsTreatA).
#	suffix: an optional prefix (such as the treatment, eg WTvsTreatA).
#	fdr: add the fdr in parentheses to the end? eg " (0.001)"
#	rank: add the rank at the start? eg "12: "
#	direction: add UP or DOWN to the end?
#	maxlen: the max number of characters in the geneset name (before we add rank & fdr & direction to it)
#
# Value:
#	a GSEA list with modified names.
#
# Mark Cowley, 2009-10-14
#


##' Rename the genesets from within a GSEA object.
##' 
##' Rename genesets using a prefix and suffix, such that COWLEY_UP becomes:
##' PREFIX_COWLEY_UP_SUFFIX
##' Optionally, the rank, FDR and direction can also be included so that
##' COWLEY_UP becomes:
##' 1: PREFIX_COWLEY_UP_SUFFIX UP (0.002)
##' 
##' @param x a GSEA list, optionally with an edb. see import.gsea(...,
##'   edb=TRUE)
##' @param prefix an optional prefix (such as the treatment, eg WTvsTreatA).
##' @param suffix an optional prefix (such as the treatment, eg WTvsTreatA).
##' @param fdr add the fdr in parentheses to the end? eg " (0.001)"
##' @param rank add the rank at the start? eg "12: "
##' @param direction add UP or DOWN to the end?
##' @param maxlen the max number of characters in the geneset name (before we
##'   add rank & fdr & direction to it)
##' @return a GSEA list with modified names.
##' @author Mark Cowley, 2009-10-14
##' @export
gsea.rename.genesets <- function(x, prefix=NULL, suffix=NULL, fdr=FALSE, rank=FALSE, direction=FALSE, maxlen=60) {

	old.names <- x$tt$NAME
	new.names <- .gsea.rename.genesets(old.names, x$tt, direction=direction, fdr=fdr, prefix=prefix, suffix=suffix, rank=rank, maxlen=maxlen)

	x$tt$NAME <- new.names[match(old.names, x$tt$NAME)]
	names(x$leading.edge) <- new.names[match(old.names, names(x$leading.edge))]
	if( "edb" %in% names(x) ) {
		# x$edb <- gsea.rename.genesets.edb(x$edb, prefix=prefix, suffix=suffix, fdr=fdr, rank=rank, direction=direction, maxlen=maxlen)
		x$edb <- gsea.rename.genesets.edb(x$edb, old.names, new.names)
		names(x$gmt) <- new.names[match(old.names, names(x$gmt))]
	}

	x
}


# Internal function used to rename geneset names, based on a prefix and suffix, and optionally, their rank, qvalue, and direction.
#
# Parameters:
#	genesets: a vector of geneset names
#	tt: a gsea top table. see import.gsea.topTable
#	prefix: an optional prefix (such as the treatment, eg WTvsTreatA). Separated by a single space.
#	suffix: an optional prefix (such as the treatment, eg WTvsTreatA). Separated by a single space.
#	fdr: add the fdr in parentheses to the end? eg " (0.001)"
#	rank: add the rank at the start? eg "12: "
#	direction: add UP or DOWN to the end?
#	maxlen: the max number of characters in the geneset name (before we add rank & fdr & direction to it)
#
# Value:
#	a vector of geneset names, same length and order as input, something like this:
#	  "1: HSA04512_ECM_RECEPTOR_INTERACTION (0.00952) up"
#     "3: HSA04060_CYTOKINE_CYTOKINE_RECEPTOR_I... (0.00645) up"
#
# Mark Cowley, 2009-09-03
.gsea.rename.genesets <- function(genesets, tt, prefix=NULL, suffix=NULL, fdr=FALSE, rank=FALSE, direction=FALSE, maxlen=60) {
	idx <- match(genesets, tt$NAME)
	if( any(is.na(idx)) )
		stop("Not all genesets found in the tt$NAME column.\n")

	tt <- tt[idx, c("RANK", "NAME", "FDR.q.val", "DIRECTION")]
	
	if( any(nchar(tt$NAME) > maxlen) ) {
		idx <- nchar(tt$NAME) > maxlen
		tt$NAME[idx] <- paste(str.left(tt$NAME[idx], maxlen-3), "...", sep="")
	}
	
	# V1: not flexible enough.
	# res <- apply(tt[,1:3], 1, function(x) {
	# 	sprintf("%s: %s (%.3f)", as.character(x[1]), x[2], as.numeric(x[3]))
	# })
	# if( direction )
	# 	res <- paste(res, tt$DIRECTION)	

	# V2
	res <- tt$NAME
	if( !is.null(prefix) )
		res <- sprintf("%s_%s", prefix, res)
	if( !is.null(suffix) )
		res <- sprintf("%s_%s", res, suffix)
	if( fdr )
		res <- sprintf("%s (%.3f)", res, as.numeric(tt$FDR.q.val))
	if( rank )
		res <- sprintf("%d: %s", as.numeric(tt$RANK), res)
	if( direction )
		res <- sprintf("%s %s", res, toupper(tt$DIRECTION))
	
	res
}



# Parse an edb XML object, and rename the genesets with optional prefix and suffix.
# NB: any spaces in pre/suffix will be replaced with "_".
# I suggest separating words with "_". untested using :, ;, -, '.'
#
# Parameters:
#	edb: an XMLDocument. See import.gsea.edb(..., edb=TRUE)
#	old.names: a vector of original geneset names. these should be found within the edb names
#	new.names: the new names, in the same order as the old names.
#
# Details:
#	This is intended to be used by gsea.rename.genesets, which should calculate old.names and new.names. Entries in edb, given that they are xml are in no strict order (although they seem to be ordered from most -ve NES to most +ve NES). Thus you need to supply the old.names is so that we can ensure the entries in edb get renamed correctly.
#
# Value:
#	An mutable XMLInternalDOM (ie an XML tree much like an XMLDocument but which can be modified) which is a copy of the original edb XMLDocument.
#
# Mark Cowley, 2009-10-06
#


##' Parse an edb XML object, and rename the genesets with optional prefix and
##' suffix.
##' 
##' This is intended to be used by gsea.rename.genesets, which should calculate
##' old.names and new.names. Entries in edb, given that they are xml are in no
##' strict order (although they seem to be ordered from most -ve NES to most
##' +ve NES). Thus you need to supply the old.names is so that we can ensure
##' the entries in edb get renamed correctly.
##' 
##' @param edb an XMLDocument. See import.gsea.edb(..., edb=TRUE)
##' @param old.names a vector of original geneset names. these should be found
##'   within the edb names
##' @param new.names the new names, in the same order as the old names.
##' @return An mutable XMLInternalDOM (ie an XML tree much like an XMLDocument
##'   but which can be modified) which is a copy of the original edb
##'   XMLDocument.
##' @author Mark Cowley, 2009-10-06
##' @export
gsea.rename.genesets.edb <- function(edb, old.names, new.names) {

	r <- xmlRoot(edb, skip=TRUE)

	edb.names <- .edb.names(edb)
	geneset.prefix <- rep("gene_sets.gmt", length(edb.names))
	
	# # need to ensure that the names in the edb are in the same order as the old.names and thus new.names.
	# edb.names <- as.character(xmlApply(r, function(node) { xmlAttrs(node)["GENESET"] }))
	# tmp <- strsplit(edb.names, "#")
	# geneset.prefix <- sapply(tmp, "[", 1) # the LHS
	# edb.names <- sapply(tmp, "[", 2) # the RHS

	idx <- match(edb.names, old.names)
	if(any(is.na(idx)))
		stop("Some of the geneset names within edb are not found in your old.names vector.")

	old.names <- old.names[idx]
	new.names <- new.names[idx]
	new.names <- paste(geneset.prefix, new.names, sep="#")

	suppressWarnings(res <- xmlTree("EDB", attrs=xmlAttrs(r))) # produces these warnings: 
	# In xmlRoot.XMLInternalDocument(currentNodes[[1]]) : empty XML document
	# even in the example(xmlTree)

	for(i in 1:xmlSize(r)) {
		a <- xmlAttrs(r[[i]])
		a["GENESET"] <- new.names[i]

		res$addNode("DTG", attrs=a)
	}

	res
}


# Private function to get the edb GENESET names from within an edb object.
#
#
.edb.names <- function(edb) {
	r <- xmlRoot(edb, skip=TRUE)

	# need to ensure that the names in the edb are in the same order as the old.names and thus new.names.
	edb.names <- as.character(xmlApply(r, function(node) { xmlAttrs(node)["GENESET"] }))
	tmp <- strsplit(edb.names, "#")
	# geneset.prefix <- sapply(tmp, "[", 1) # the LHS
	edb.names <- sapply(tmp, "[", 2) # the RHS
	
	edb.names
}

# # Parse an edb XML object, and rename the genesets with optional prefix and suffix.
# # NB: any spaces in pre/suffix will be replaced with "_".
# # I suggest separating words with "_". untested using :, ;, -, '.'
# #
# # Parameters:
# #	edb: an XMLDocument. See import.gsea.edb
# #	prefix: A character(1). eg "PRE_", or "SET_A_"
# #	suffix: A character(1). eg "_POST", or "_SET_B"
# #	fdr: add the fdr in parentheses to the end? eg " (0.001)"
# #	rank: add the rank at the start? eg "12: "
# #	direction: add UP or DOWN to the end?
# #	maxlen: the max number of characters in the geneset name (before we add rank & fdr & direction to it)
# #
# # Details:
# #	The entries in edb files are unaware of their ranking (which can be determined from the NES score.)
# #
# # Value:
# #	An mutable XMLInternalDOM (ie an XML tree much like an XMLDocument but which can be modified) which is a copy of the original edb XMLDocument.
# #
# # Mark Cowley, 2009-10-06
# #
# gsea.rename.genesets.edb <- function(edb, prefix=NULL, suffix=NULL, fdr=FALSE, rank=FALSE, direction=FALSE, maxlen=60) {
# 
# 	r <- xmlRoot(edb, skip=TRUE)
# 
# 	#
# 	# I want to use the .gsea.rename.genesets to change my geneset names.
# 	# This will be FAR better for consistency.
# 	# Thus I need to make a tt object, containing these columns at a minimum:
# 	# c("RANK", "NAME", "FDR.q.val", "DIRECTION")
# 	#
# 	
# 	elements <- xmlApply(r, function(node) { xmlAttrs(node)[c("GENESET", "NES", "FDR")] })
# 	old.names <- as.character(sapply(elements, "[", 1))
# 	tmp <- strsplit(old.names, "#")
# 	geneset.prefix <- sapply(tmp, "[", 1) # the LHS
# 	old.names <- sapply(tmp, "[", 2) # the RHS
# 
# 	nes.scores <- as.numeric(sapply(elements, "[", 2))
# 	fdr.vals <- as.numeric(sapply(elements, "[", 3))
# 
# 	# edb entries don't know thir ranking relative to each other.
# 	# Determine the ranking based on the NES scores of all entries....
# 	ranks <- rep(NA, length(nes.scores))
# 	ranks[nes.scores>0] <- rank(1/nes.scores[nes.scores>0])
# 	ranks[nes.scores<0] <- rank(-1/nes.scores[nes.scores<0])
# 
# 	direc <- rep("UP", length(nes.scores))
# 	direc[nes.scores<0] <- "DOWN"
# 	
# 	# make the tt object...
# 	tt <- data.frame(RANK=ranks, NAME=old.names, NES=nes.scores, FDR.q.val=fdr.vals, DIRECTION=direc, stringsAsFactors=FALSE)
# 	
# 	# determine the new names for the genesets.
# 	new.names <- .gsea.rename.genesets(genesets=old.names, tt=tt, prefix=prefix, suffix=suffix, fdr=fdr, rank=rank, direction=direction, maxlen=maxlen)
# 	new.names <- paste(geneset.prefix, new.names, sep="#")
# 
# 	suppressWarnings(res <- xmlTree("EDB", attrs=xmlAttrs(r))) # produces these warnings: 
# 	# In xmlRoot.XMLInternalDocument(currentNodes[[1]]) : empty XML document
# 	# even in the example(xmlTree)
# 
# 	for(i in 1:xmlSize(r)) {
# 		a <- xmlAttrs(r[[i]])
# 		a["GENESET"] <- new.names[i]
# 		# geneset.path <- strsplit(a[3], "#")[[1]]
# 		# geneset <- geneset.path[2]
# 		# if( !is.null(prefix) )
# 		# 	geneset <- sprintf("%s_%s", prefix, geneset)
# 		# if( !is.null(suffix) )
# 		# 	geneset <- sprintf("%s_%s", geneset, suffix)
# 		# if( fdr )
# 		# 	geneset <- sprintf("%s (%.3f)", geneset, as.numeric(a["FDR"]))
# 		# if( rank ) { 
# 		# 	geneset <- sprintf("%d: %s", ranks[i], geneset)
# 		# }
# 		# if( direction ) {
# 		# 	direc <- ifelse(as.numeric(a["NES"])>0,"UP","DOWN")
# 		# 	geneset <- sprintf("%s %s", geneset, direc)
# 		# }
# 		# geneset.path[2] <- geneset
# 		# geneset.path <- paste(geneset.path, collapse="#")
# 		# a[3] <- geneset.path
# 
# 		res$addNode("DTG", attrs=a)
# 	}
# 
# 	res
# }
