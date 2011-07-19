#######   WARNING  ############
#
# This is probably defunct, thanks to gsea.merge, which merges GSEA objects that are already loaded into R.
# This code is less well tested than the gsea.merge
#
# Mark Cowley, 2009-10-16
#




#' Merge 2 gsea runs, by referencing only the GSEA directories.
#' 
#' This only merges the contents of the edb directory. Merging the individual
#' html geneset reports is very complex...
#' Why?
#' so you can examine them both together in the Leading Edge analysis tool of
#' course!
#' See also gsea.merge
#' 
#' @param gsea.dir1 paths to two GSEA directories (the ones that contain
#'   index.html)
#' @param gsea.dir2 paths to two GSEA directories (the ones that contain
#'   index.html)
#' @param prefix both should be character(2) indicating how you want the
#'   genesets to be renamed so that they reflect which experimental run they
#'   came from.
#' @param suffix both should be character(2) indicating how you want the
#'   genesets to be renamed so that they reflect which experimental run they
#'   came from.
#' @param outdir The path to an edb directory that will be created.
#' @author Mark Cowley, 2009-10-06
#' @export
gsea.merge.dirs <- function(gsea.dir1, gsea.dir2, prefix, suffix, outdir) {
	stopifnot( !missing(outdir) )
	stopifnot( is.gsea.dir(gsea.dir1) )
	stopifnot( is.gsea.dir(gsea.dir2) )
	if( missing(prefix) ) prefix <- c("", "")
	if( missing(suffix) ) suffix <- c("", "")

	if( !file.exists(outdir) ) dir.create( outdir )
	outdir <- file.path(outdir, "edb")
	if( !file.exists(outdir) ) dir.create( outdir )

	
	# copy the rnk file across
	cat("Copying rnk file...")
	rnk1 <- dir(file.path(gsea.dir1, "edb"), pattern=".*rnk", full=TRUE)
	file.copy(rnk1, outdir, overwrite=TRUE)
	cat("\n")
	
	# merge and rename the gmt file "gene_sets.gmt"
	cat("Merging gmt files..")
	gmt1 <- import.gsea.gmt(file.path(gsea.dir1, "edb", "gene_sets.gmt"))
	gmt2 <- import.gsea.gmt(file.path(gsea.dir2, "edb", "gene_sets.gmt"))
	names(gmt1) <- paste(prefix[1], names(gmt1), suffix[1], sep="")
	names(gmt2) <- paste(prefix[2], names(gmt2), suffix[2], sep="")
	gmt <- lbind(gmt1, gmt2)
	export.gsea.gmt(gmt, file.path(outdir, "gene_sets.gmt"))
	cat("\n")
	
	# rename, merge and export the edb files
	cat("Merging edb files...")
	edb <- gsea.merge.edb.dirs(gsea.dir1, gsea.dir2, prefix, suffix)
	export.gsea.edb(edb, file.path(outdir, "results.edb"))
	cat("\n")
}


#' Function to merge the edb files from 2 GSEA runs.
#' 
#' Why?
#' so you can examine them both together in the Leading Edge analysis tool of
#' course!
#' Warnings:
#' The resultant edb must have the same RANKED_LIST attribute for all nodes.
#' Thus the RANKED_LIST
#' attribute for the 2nd GSEA run is changed to the name of the first rnk
#' which is obviously a bit
#' of a hack.
#' 
#' @param gsea.dir1 paths to two GSEA directories (the ones that contain
#'   index.html)
#' @param gsea.dir2 paths to two GSEA directories (the ones that contain
#'   index.html)
#' @param prefix both should be character(2) indicating how you want the
#'   genesets to be renamed so that they reflect which experimental run they
#'   came from.
#' @param suffix both should be character(2) indicating how you want the
#'   genesets to be renamed so that they reflect which experimental run they
#'   came from.
#' @return none. merges 2 directories
#' @author Mark Cowley, 2009-10-06
#' @export
gsea.merge.edb.dirs <- function(gsea.dir1, gsea.dir2, prefix, suffix) {
	require(XML)
	edb1 <- import.gsea.edb(gsea.dir1)
	edb2 <- import.gsea.edb(gsea.dir2)
	
	edb1 <- gsea.rename.genesets.edb(edb1, prefix[1], suffix[1])
	edb2 <- gsea.rename.genesets.edb(edb2, prefix[2], suffix[2])
	
	r1 <- xmlRoot(edb1, skip=TRUE)
	r2 <- xmlRoot(edb2, skip=TRUE)
	
	a <- NULL
	suppressWarnings(res <- xmlTree("EDB", attrs=xmlAttrs( r1 )))
	for(i in 1:xmlSize(r1)) {
		a <- xmlAttrs(r1[[i]])
		res$addNode("DTG", attrs=a)
	}
	rnk <- a[1]
	for(i in 1:xmlSize(r2)) {
		a <- xmlAttrs(r2[[i]])
		a[1] <- rnk # GSEA complains if you have different rnk files present.
		res$addNode("DTG", attrs=a)
	}

	res 
}
