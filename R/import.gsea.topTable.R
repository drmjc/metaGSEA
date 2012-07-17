#' Import the tabular results from a GSEA PreRanked run
#' 
#' Import the tabular results from a GSEA PreRanked run - ie the GSEA top
#' table.
#' 
#' Takes the pos and neg results, and merges them into a single table, adding
#' a DIRECTION column.
#' Then re-sort the rows by NES (default), absolute NES, FDR, P or SIZE
#' (largest->smallest)
#' 
#' @param dir the dir that contains the index.html
#' @param sort.by which column to sort by.
#' @param NES decreasing NES score
#' @param absNES decreasing absolute NES score
#' @param FDR increasing FDR, breaking ties by largest absolute NES
#' @param P increasing nomimal P value, breaking ties by largest absolute NES
#' @param SIZE decreasing gene set size
#' 
#' @return a single \code{data.frame} with all genesets from within a collection.
#'   Unlike the raw GSEA output, this contains a \code{DIRECTION} column (up/down)
#'   and a \code{LEADING.EDGE.SIZE} column which is the number of genes in the
#'   leading edge +/- 1 due to rounding errors (this data comes from the
#'   \code{tags} field in the \code{LEADING.EDGE} column), and a \code{RANK.IN.REPORT} column
#'   which shows the rank in the original list
#' 
#' @author Mark Cowley, 2009-04-28
#' @export
import.gsea.topTable <- function(dir, sort.by=c("NES", "absNES", "FDR", "P", "SIZE")) {
	if( !file.exists(file.path(dir, "index.html")) ) {
		stop("Must specify the top-level dir that contains index.html.\n")
	}

	####################
	# determine the 2 class names that have been compared.
	# for GseaPreranked, this will be na_pos and na_neg
	# for Gsea, we need to look inside the index.html file to find out.
	#
	rpt <- import.gsea.rpt(dir)
	mode <- rpt$producer_class
	if( mode == "xtools.gsea.GseaPreranked") {
		classes <- c("na_pos", "na_neg")
	}
	else if( mode == "xtools.gsea.Gsea" ) { # look inside the index.html to find the order.
		classes <- .gsea.get.classes.index.html(file.path(dir, "index.html"))
	}
	patterns <- sprintf("gsea_report_for_%s.*\\.xls", classes)
	files <- c( dir(dir, pattern=patterns[1], full=TRUE), dir(dir, pattern=patterns[2], full.names=TRUE) )
	if( !all(file.exists(files)) ) {
		stop("Can't find the 2 top table xls files (gsea_report_for_.*.xls).\n")
	}
	# done.
	#####################
	
	####################
	# import the 2 files, and then combine & re-sort.
	#
	tt0 <- .import.gseaPreranked.topTable.file(files[1])
	tt1 <- .import.gseaPreranked.topTable.file(files[2])
	tt0$DIRECTION <- "up"
	tt0$ENRICHED.PHENOTYPE <- classes[1]
	tt1$DIRECTION <- "down"
	tt1$ENRICHED.PHENOTYPE <- classes[2]
	
	tt <- rbind(tt0, tt1)
	tt <- move.column(tt, c("DIRECTION", "ENRICHED.PHENOTYPE"), "ES")
	
	# re-sort the table
	tt <- gsea.topTable.sort(tt, sort.by=sort.by[1])
	#####################

	return( tt )
	
}

gsea.topTable.sort <- function(tt, sort.by=c("NES", "absNES", "FDR", "P", "SIZE")) {
	sort.by <- sort.by[1]
	if( sort.by == "NES" )
		tt <- tt[order(tt$NES, decreasing=TRUE), ]
	else if( sort.by == "absNES" )
		tt <- tt[order(abs(tt$NES), decreasing=TRUE), ]
	else if( sort.by == "FDR" )
		tt <- tt[order(tt$FDR.q.val, -abs(tt$NES), decreasing=FALSE), ]
	else if( sort.by == "P" )
		tt <- tt[order(tt$NOM.p.val, -abs(tt$NES), decreasing=FALSE), ]
	else if( sort.by == "SIZE" )
		tt <- tt[order(tt$SIZE, decreasing=TRUE), ]
	rownames(tt) <- 1:nrow(tt)

	return( tt )
	
}

# When running GSEA, you don't know which is class0 or class 1 unless you look inside the index.html.
# This function returns the classes in the relevant order.
#
# Mark Cowley, 2009-12-10
.gsea.get.classes.index.html <- function(f) {
	tmp <- suppressWarnings(readLines(f))
	tmp <- tmp[grep("Enrichment in phenotype", tmp)]
	class0 <- sub(" \\(.*", "", sub("</div><div><h4>Enrichment in phenotype: <b>", "", tmp))
	class1 <- sub(" \\(.*", "", sub("^</div><div><h4>Enrichment.+<h4>Enrichment in phenotype: <b>", "", tmp))
	res <- c(class0, class1)
	return(res)
}

# Private function to import a pos/neg summary table (in xls)
#
# Mark Cowley, 2009-12-10
.import.gseaPreranked.topTable.file <- function(f) {
	res <- read.delim(f)
	if( nrow(res) == 0 )
		return( NULL )

	#
	# remove junk columns, add DIRECTION and add LEADING.EDGE.SIZE columns
	#
	cols <- setdiff(colnames(res), c("GS.br..follow.link.to.MSigDB", "GS.DETAILS", "X"))
	res <- res[,cols]
	tags <- res$LEADING.EDGE
	tags <- sub("%.*", "", tags)
	tags <- sub("tags=", "", tags)
	tags <- as.numeric(tags) / 100
	res$LEADING.EDGE.SIZE <- round(res$SIZE * tags)
	################
	# sometimes you can get a NES of '' and NOM.p.val of '' which gives you an NA.
	res$ES[is.na(res$ES)] <- 0.0
	res$NES[is.na(res$NES)] <- 0.0
	res$NOM.p.val[is.na(res$NOM.p.val)] <- 1.0
	res$FDR.q.val[is.na(res$FDR.q.val)] <- 1.0
	res$FWER.p.val[is.na(res$FWER.p.val)] <- 1.0
	# res <- res[order(res$NES, decreasing=(i==1)), ]
	################
	
	res$RANK <- 1:nrow(res) # make sure rank is in the same order as it was in the html table. ie unchanged.
	res <- move.column(res, "RANK", "ES")
	
	return( res )
}



# .import.gseaGSEA.topTable <- function(dir, sort.by=c("NES", "absNES", "FDR", "P", "SIZE")) {
# 	f <- dir(dir, pattern="gsea_report_for_.*\\.xls", full.names=TRUE)
# 	classes <- sub("_.*", "", sub("gsea_report_for_", "", f))
# 	
# }
# 
# .import.gseaPreranked.topTable <- function(dir, sort.by=c("NES", "absNES", "FDR", "P", "SIZE")) {
# 	tt.list <- .import.gseaPreranked.topTable.posneg(dir)
# 	for(i in 1:2) {
# 		try(tt.list[[i]]$DIRECTION <- c("up", "down")[i], silent=TRUE) # ignore the error if tt.list[[i]] has 0 rows
# 	}
# 	tt <- rbind(tt.list[[1]], tt.list[[2]])
# 	
# 	# re-sort the table
# 	sort.by <- sort.by[1]
# 	if( sort.by == "NES" )
# 		tt <- tt[order(tt$NES, decreasing=TRUE), ]
# 	else if( sort.by == "absNES" )
# 		tt <- tt[order(abs(tt$NES), decreasing=TRUE), ]
# 	else if( sort.by == "FDR" )
# 		tt <- tt[order(tt$FDR.q.val, -abs(tt$NES), decreasing=FALSE), ]
# 	else if( sort.by == "P" )
# 		tt <- tt[order(tt$NOM.p.val, -abs(tt$NES), decreasing=FALSE), ]
# 	else if( sort.by == "SIZE" )
# 		tt <- tt[order(tt$SIZE, decreasing=TRUE), ]
# 	rownames(tt) <- 1:nrow(tt)
# 
# 	tt <- move.column(tt, "DIRECTION", "ES")
# 	# done
# 	return( tt )
# }
# 
# #
# # Private function to import GSEA results from a single run into a list with pos/neg elements.
# #
# # Import the pos and neg xls files from a single GSEA rnk vs one gmt file.
# # Parameters:
# #	dir: the path to the top level directory, eg "./GSEA/c1_all.GseaPreranked.1220343365997"
# #
# # Value:
# #	a list with 2 data.frames called pos and neg, for the up- and down-regulated results
# #
# # Mark Cowley, 2008-09-03
# # import.gsea.1 <- function(dir) {
# .import.gseaPreranked.topTable.posneg <- function(dir) {
# 	res <- list()
# 	f <- dir(dir, pattern="gsea_report.*_pos_.*xls", full.names=TRUE)
# 	if( length(f) == 1 && file.exists(f) ) {
# 		res$pos <- read.delim(f)
# 	}
# 	else {
# 		res$pos <- NA
# 	}
# 	
# 	f <- dir(dir, pattern="gsea_report.*_neg_.*xls", full.names=TRUE)
# 	if( length(f) == 1 && file.exists(f) ) {
# 		res$neg <- read.delim(f)
# 	}
# 	else {
# 		res$neg <- NA
# 	}
# 
# 	#
# 	# remove junk columns, add DIRECTION and add LEADING.EDGE.SIZE columns
# 	#
# 	cols <- setdiff(colnames(res$pos), c("GS.br..follow.link.to.MSigDB", "GS.DETAILS", "X"))
# 	for(i in 1:2) {
# 		if( nrow(res[[i]]) == 0 )
# 			next
# 		res[[i]] <- res[[i]][,cols]
# 		tags <- res[[i]]$LEADING.EDGE
# 		tags <- sub("%.*", "", tags)
# 		tags <- sub("tags=", "", tags)
# 		# tags <- str.right(str.left(res[[i]]$LEADING.EDGE, 7), 2)
# 		tags <- as.numeric(tags) / 100
# 		res[[i]]$LEADING.EDGE.SIZE <- round(res[[i]]$SIZE * tags)
# 		################
# 		# sometimes you can get a NES of '' and NOM.p.val of '' which gives you an NA.
# 		res[[i]]$ES[is.na(res[[i]]$ES)] <- 0.0
# 		res[[i]]$NES[is.na(res[[i]]$NES)] <- 0.0
# 		res[[i]]$NOM.p.val[is.na(res[[i]]$NOM.p.val)] <- 1.0
# 		res[[i]]$FDR.q.val[is.na(res[[i]]$FDR.q.val)] <- 1.0
# 		res[[i]]$FWER.p.val[is.na(res[[i]]$FWER.p.val)] <- 1.0
# 		# res[[i]] <- res[[i]][order(res[[i]]$NES, decreasing=(i==1)), ]
# 		################
# 		
# 		res[[i]]$RANK <- 1:nrow(res[[i]]) # make sure rank is in the same order as it was in the html table. ie unchanged.
# 		res[[i]] <- move.column(res[[i]], "RANK", "ES")
# 	}
# 	
# 	return( res )
# }


