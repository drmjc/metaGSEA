#' Export a gsea object to the filesystem.
#' 
#' This is mostly used to export the edb portion of a gsea run, since copying
#' over the individual geneset reports etc sounds far too tedious, plus these
#' elements are not stored inside the GSEA object in R.
#' 
#' @param x a GSEA object.
#' @param location the dir that will be created. eg: \dQuote{./GSEA/c2_all}, 
#'   \dQuote{c2_all.zip}, and set \code{as.zip=TRUE}
#' @param as.zip logical: create a zip file? If \code{TRUE}, then make sure the
#'   \code{location} parameter is set to a zip file name, rather than a directory
#'   path.
#' @param overwrite logical: overwrite an existing file or directory?
#' @return none. Creates either a directory, or a zip file of a directory. 
#' It will only contain the information that was loaded by import.gsea, which means
#' that individual geneset html resports and pictures are not exported.
#' @author Mark Cowley, 2009-10-15
#' @export
export.gsea <- function(x, location, as.zip=FALSE, overwrite=FALSE) {
	if( as.zip ) {
		dir <- tempdir()
		export.gsea(x, dir, as.zip=FALSE)
		zip(dir, location)
		unlink(dir)
	}
	else {
		if( file.exists(location) ) {
			if( overwrite )
				dir.remove(location)
			else
				stop("location already exists. Either change loation, or specify overwrite=TRUE")
		}
		dir.create(location)

		write.delim(x$tt, file.path(location, paste("gsea_report_combined_", x$rpt$producer_timestamp,".xls", sep="")), row.names=FALSE, col.names=TRUE)

		if( "edb" %in% names(x)) {
			dir.create(file.path(location, "edb"))
			if( "rnk1" %in% names(x) ) {
				# then this is a merged GSEA object. see gsea.merge
				export.gsea.rnk(x$rnk1, names(x$rnk1), file.path(location, "edb", basename(x$rpt1$collapsed_rnk)), resort=FALSE)
				export.gsea.rnk(x$rnk2, names(x$rnk2), file.path(location, "edb", basename(x$rpt2$collapsed_rnk)), resort=FALSE)
			}
			else {
				export.gsea.rnk(x$rnk, names(x$rnk), file.path(location, "edb", basename(x$rpt$collapsed_rnk)), resort=FALSE)
			}
			export.gsea.gmt(x$gmt, file.path(location, "edb", "gene_sets.gmt"))
			export.gsea.edb(x$edb, file.path(location, "edb", "results.edb"))
		}
	}
}
