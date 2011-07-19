# Import the output produced by gsea.compare.runs.sh
#
# Parameters:
#	xls: an XLS file produced by running import.gsea.compare.runs.sh, or
# 		 import.gsea.compare.runs.1gmt.sh
#	gmt: either an integer indicating which sheet of the xls to import, or
#		the name of the gmt file. use the latter if gsea was run on all 15 of the
#		gmt categories listed in .GSEA.GMT.CATEGORIES
#
# Examples:
#	# not run
# 	# gseacmp <- import.gsea.compare.runs("GSEA/gsea.compare.runs.xls")
# 	# gseacmp <- import.gsea.compare.runs("GSEA/gsea.compare.runs.xls", c("c2_all"))
# 	# gseacmp <- import.gsea.compare.runs("GSEA/gsea.compare.runs.xls", c("c2_all", "c5_bp"))
# 	# gseacmp <- import.gsea.compare.runs("GSEA/gsea.compare.runs.xls", grep("all", .GSEA.GMT.CATEGORIES, value=TRUE))
# 
# Mark Cowley, 2009-03-23


##' Import the output produced by gsea.compare.runs.sh
##' 
##' @param xls an XLS file produced by running import.gsea.compare.runs.sh, or
##'   import.gsea.compare.runs.1gmt.sh
##' @param gmt either an integer indicating which sheet of the xls to import,
##'   or the name of the gmt file. use the latter if gsea was run on all 15 of
##'   the gmt categories listed in .GSEA.GMT.CATEGORIES
##' @author Mark Cowley, 2009-03-23
##' @examples
##' # not run
##' # gseacmp <- import.gsea.compare.runs("GSEA/gsea.compare.runs.xls")
##' # gseacmp <- import.gsea.compare.runs("GSEA/gsea.compare.runs.xls", c("c2_all"))
##' # gseacmp <- import.gsea.compare.runs("GSEA/gsea.compare.runs.xls", c("c2_all", "c5_bp"))
##' # gseacmp <- import.gsea.compare.runs("GSEA/gsea.compare.runs.xls", grep("all", .GSEA.GMT.CATEGORIES, value=TRUE))
##' @export
import.gsea.compare.runs <- function(xls, gmt=.GSEA.GMT.CATEGORIES) {
	if( is.null(gmt) )
		gmt <- .GSEA.GMT.CATEGORIES

	if( length(gmt) > 1 ) {
		res <- list()
		for(gmt in .GSEA.GMT.CATEGORIES) {
			res[[gmt]] <- import.gsea.compare.runs(xls, gmt)
		}
		return( res )
	}
	else {
		if( is.character(gmt) )
			gmt <- match(gmt, .GSEA.GMT.CATEGORIES)
		res <- read.xls(xls, sheet=gmt)
		return( res )
	}
}
