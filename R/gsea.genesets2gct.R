#
# Create gct files for all gseaPreranked runs comparing 1 rnk to lots of gmt's.
#
# By default, the leading edge genes are identified by the 'leading.edge' slot from
# running import.gsea(); this then creates an edb/leading_edge.gmt file.
# From a gct.file at the unique genesymbol level (see gsea.gct.probes2genes), the
# relevant rows are retrieved, keeping the same order as the genes in the leading 
# edge itself.
#
# Parameters:
#	gsea.dir: either a single gsea directory, or a parent-dir that contains lots 
# 		of gsea dirs, (each vs a different collection of genesets)
#	gct.file: a unique gene-level gct file (very important to see gsea.gct.probes2genes)
#	leading.edge: if TRUE then only the leading edge genes are exported,
#		or else, the genes in the same order as given by edb/gene_sets.gmt are given.
#	which.sets: do you want to make a gct for every single leading edge ("all"), or just those
#		that have been created in the dir ("best") which are usually the top 20 or 50.
#
# Mark Cowley, 2009-07-27


##' Create gct files for all gseaPreranked runs comparing 1 rnk to lots of
##' gmt's.
##' 
##' By default, the leading edge genes are identified by the 'leading.edge'
##' slot from
##' running import.gsea(); this then creates an edb/leading_edge.gmt file.
##' From a gct.file at the unique genesymbol level (see gsea.gct.probes2genes),
##' the
##' relevant rows are retrieved, keeping the same order as the genes in the
##' leading
##' edge itself.
##' 
##' @param gsea.dir either a single gsea directory, or a parent-dir that
##'   contains lots of gsea dirs, (each vs a different collection of genesets)
##' @param gct.file a unique gene-level gct file (very important to see
##'   gsea.gct.probes2genes)
##' @param leading.edge if TRUE then only the leading edge genes are exported,
##'   or else, the genes in the same order as given by edb/gene_sets.gmt are
##'   given.
##' @param which.sets do you want to make a gct for every single leading edge
##'   ("all"), or just those that have been created in the dir ("best") which
##'   are usually the top 20 or 50.
##' @author Mark Cowley, 2009-07-27
##' @export
gsea.genesets2gct <- function(gsea.dir, gct.file, leading.edge=TRUE, which.sets=c("best", "all")) {
	if( !is.gsea.dir(gsea.dir) && is.dir(gsea.dir) ) {
		subdirs <- dir(gsea.dir, full=TRUE)
		if( any(is.gsea.dir(subdirs)) ) {
			subdirs <- subdirs[is.gsea.dir(subdirs)]
			cat(sprintf("Found %d GSEA subdirs.\n", length(subdirs)))
			for(subdir in subdirs) {
				cat("Generating gct's for: ", basename(subdir), "...")
				gsea.genesets2gct(gsea.dir=subdir, gct.file=gct.file, leading.edge=leading.edge, which.sets=which.sets)
				cat("\n")
			}
		}
		return(TRUE)
	}
	stopifnot( is.gsea.dir(gsea.dir) )
	
	gsea <- import.gsea(gsea.dir)
	gct <- import.gsea.gct(gct.file, check.names=FALSE)
	rownames(gct) <- toupper(rownames(gct))

	if( leading.edge ) {
		gmt <- gsea$leading.edge
		f <- file.path(gsea.dir, "edb", "leading_edge.gmt")
		export.gsea.gmt(gsea$leading.edge, f)
	}
	else {
		f <- file.path(gsea.dir, "edb", "gene_sets.gmt")
		gmt <- import.gsea.gmt(f)
		# reorder the genes to match the order in the rnk file.
		rnk <- import.gsea.rnk(gsea.dir)
		gmt <- gsea.reorder.gmt.by.rnk(gmt, rnk)
	}

	gct.dir <- gsea.dir

	which.sets <- which.sets[1]
	if( which.sets == "all" )
		genesets <- names(gmt)
	else if( which.sets == "best" ) {
		# make one gct file per xls/html file that's been made.
		xls.files <- dir(gsea.dir, pattern=".*xls")
		avoid <- c( "gene_set_sizes.xls",
					grep("gsea_report", xls.files, value=TRUE),
					grep("ranked_gene_list", xls.files, value=TRUE))
		xls.files <- setdiff(xls.files, avoid)
		genesets <- sub(".xls", "", xls.files)
	}
	
	# gct.dir <- file.path(gsea.dir, "LeadingEdge.gct")
	# if( !file.exists(gct.dir) )
	# 	dir.create(gct.dir)
	
	for(i in 1:length(genesets)) {
		gset <- genesets[i]
		genes <- gmt[[gset]]
		idx <- match(genes, rownames(gct))
		if( any(is.na(idx)) ){
			cat("The following genes(s) were missing from ", gset, ": ", tocsv(genes[is.na(idx)]), ".\n")
			idx <- na.rm(idx)
		}
		description <- gct$Description[idx]
		data <- gct[idx, 3:ncol(gct)]
		f <- file.path(gct.dir, paste(gset, ".gct", sep=""))
		export.gsea.gct(data, description=description, file=f)
	}
}
