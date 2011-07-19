# Import the rpt file from a GSEA, or GseaPreRanked analysis run
#
# If the GSEA was run on a GenePattern server, then some of the values (eg chip, gmx, cls/rnk) will point to
# locations that are unavailable on the end users machine. chip and gmx files are replaces with URL's to Broad Institute. cls, rnk, out, file references are changed to locations relative to their current location on the end users hard drive.
#
# Parameters:
#	x: the path to either a rpt file, or the top-level dir that contains index.html.
#
# Mark Cowley, 2009-04-28
# 2009-12-11: major mods to allow GSEA and GseaPreRanked rpt's to be imported that may have been run on a GenePattern server.
#


##' Import the rpt file from a GSEA, or GseaPreRanked analysis run
##' 
##' If the GSEA was run on a GenePattern server, then some of the values (eg
##' chip, gmx, cls/rnk) will point to
##' locations that are unavailable on the end users machine. chip and gmx files
##' are replaces with URL's to Broad Institute. cls, rnk, out, file references
##' are changed to locations relative to their current location on the end
##' users hard drive.
##' 
##' @param x the path to either a rpt file, or the top-level dir that contains
##'   index.html.
##' @author Mark Cowley, 2009-04-28
##' @export
import.gsea.rpt <- function(x) {
	if( is.gsea.dir(x) ) {
		f <- dir(x, pattern=".*\\.rpt$", full=TRUE)
		return( import.gsea.rpt(f) )
	}
	else if( is.file(x) ) {
		tmp <- readLines(x)
		tmp <- tmp[nchar(tmp) > 0]
		tmp <- sub("param\t", "", tmp)
		tmp <- strsplit(tmp, "\t")
		rpt <- lapply(tmp, "[", 2)
		names(rpt) <- sapply(tmp, "[", 1)
		rpt <- rpt[order(names(rpt))]

		#########################################
		# if this has been run on a genepattern server, or the results have been moved around, 
		# then some files may not be discoverable in the dir that the data is currently sitting 
		# in. Update the links in the rpt to better represent the data state that it is currently in.
		#
		path <- path.expand(dirname(x))
		
		# if the cls file is mentioned...
		if( "cls" %in% names(rpt) ) {
			if( !file.exists(rpt$cls) && file.exists(file.path(path, "edb", rpt$cls))) {
				rpt$cls <- file.path(path, "edb", rpt$cls)
			}
			else if( file.exists(rpt$cls) ) {
				rpt$cls <- path.expand(rpt$cls)
			}
		}
		else {
			rpt$cls <- NULL
		}
		
		# if the rnk file is mentioned...
		if( "rnk" %in% names(rpt) ) {
			if ( !file.exists(rpt$rnk) && file.exists(file.path(path, "edb", rpt$rnk))) {
				rpt$rnk <- file.path(path, "edb", rpt$rnk)
			}
			else if( file.exists(rpt$rnk) ) {
				rpt$rnk <- path.expand(rpt$rnk)
			}
		}
		else {
			rpt$rnk <- NULL
		}
		
		# if the gmx is a single filename, then we need to reference that file on the Broad FTP site
		if( !file.exists(rpt$gmx) && (basename(rpt$gmx) == rpt$gmx) ) {
			url <- sprintf("http://www.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/downloads/%s", rpt$gmx)
			# if( test.broad.ftp(url) )
				rpt$gmx <- url
		}
		else if( !file.exists(rpt$gmx) ) {
			rpt$gmx <- path.expand(rpt$gmx)
		}
		
		# if the chip is a single filename, then we need to reference the file on the Broad FTP site
		if( !file.exists(rpt$chip) && (basename(rpt$chip) == rpt$chip) ) {
			url <- sprintf("ftp://gseaftp.broad.mit.edu/pub/gsea/annotations/%s", rpt$chip)
			# if( test.broad.ftp(url) )
				rpt$chip <- url
		}
		else if( file.exists(rpt$chip) ) {
			rpt$chip <- path.expand(rpt$chip)
		}

		# the outdir may be on the computation server. this should become the dir where the results are now.
		if( !file.exists(rpt$out) )
			rpt$out <- path
		else
			rpt$out <- path.expand(rpt$out)

		# the file is the index.html file which should be where the results are now.
		if( !file.exists(rpt$file) )
			rpt$file <- file.path(path, "index.html")
		else
			rpt$file <- path.expand(rpt$file)

		# rpt <- lapply(rpt, function(x) gsub("//", "/", x))
		# done
		# cat("\n")
		#
		#########################################
		
		#########################################
		# collect some extra info that may be useful for merging and exporting GSEA runs.
		#
		edb.dir <- file.path(path, "edb")
		edb <- dir(edb.dir, pattern="edb$", full=TRUE)
		collapsed_rnk <- dir(edb.dir, pattern="rnk$", full=TRUE)
		gmt <- dir(edb.dir, pattern="gmt$", full=TRUE)

		extras <- unlist(list(edb=edb, collapsed_rnk=collapsed_rnk, gmt=gmt))
		rpt <- c(rpt, extras)
		#########################################

		return( rpt )
	}
	else( stop("Must specify either the rpt file itself, or the dir containing index.html.\n") )
}
