#' Import the results from a GSEA or GseaPreRanked analysis.
#' 
#' Import the output from a single GSEA PreRanked run, including: the
#' toptable, the leading edge genes, the rnk file used, and the report which
#' summarises the analysis settings used, and optionally, the edb file.
#' 
#' A flexible function which takes a directory of GSEA results, and imports
#' the `top table', the leading edge genes, the pre-ranked file (if one was
#' used) and the report file which describes the parameters used when the GSEA
#' was run.
#' 
#' The top table is a combination of the up/down tables that GSEA makes,
#' sorted by decreasing NES, with additional `LEADING.EDGE.SIZE' and
#' `DIRECTION' columns.
#' 
#' The leading.edge is a named list, where the names are the geneset names,
#' and each list element is a vector of gene symbols, sorted in decreasing
#' significance (i.e. first gene is the most important gene for this geneset).
#' 
#' The edb file contains the machine readable output which the GSEA Desktop
#' tool can read. You should import it if you plan on subsetting or merging
#' GSEA results which are destined to be viewed in the Leading Edge Viewer
#' tool. Otherwise it is quite large, so best avoided unless you really need
#' it!
#' 
#' @param dir The path to a directory containing GSEA results. There should be
#'   an index.html file in this directory. You may specify a single directory,
#'   a vector of multiple GSEA directories, or the parent folder which
#'   contains multiple GSEA directories.
#' @param edb Import the edb file? TRUE or FALSE. See Details.
#' @return A GSEA \code{list}, with the following elements:
#'   \item{tt}{the combined (pos & neg) GSEA top table}
#'   \item{leading.edge}{a named list of leading edge genes}
#'  \item{rnk}{The pre-ranked file (if one was used)}
#' \item{rpt}{The GSEA report, listing the parameter values used during the GSEA run}
#' \item{edb}{[optional] The edb object, as an \code{XML} tree. See \code{\link{import.gsea.edb}}}
#' @author Mark Cowley
#' @seealso \code{\link{import.gsea.topTable}}
#'   \code{\link{import.gsea.leadingedge}} \code{\link{import.gsea.rpt}}
#'   \code{\link{import.gsea.rnk}} \code{\link{import.gsea.edb}}
#' @keywords IO file
#' @examples
#' # not run
#' # import.gsea("./c2_all.Gsea.1252055214188")
#' # import.gsea("./c2_all.Gsea.1252055214188", edb=TRUE)
#' # import.gsea(c("./c1_all.Gsea.1252052322484", "./c2_all.Gsea.1252055214188", "./c3_all.Gsea.1252086970993"))
#' # import.gsea("/path/to/dir/containing/lots/of/GSEA/dirs")
#' @export
import.gsea <- function(dir, edb=FALSE) {
	if( length(dir) > 1 ) {
		# then there was a vector of GSEA dirs
		dirs <- dir
		stopifnot( all(is.gsea.dir(dirs)) )
		cat(sprintf("Importing %d GSEA sub-directories.\n", length(dirs)))
		res <- list()
		for(i in 1:length(dirs)) {
			dir <- dirs[i]
			res[[i]] <- import.gsea(dir, edb=edb)
			names(res)[i] <- gsea.which.gmt(dir)
			cat(".")
		}
		if( alleq(names(res)) ) {
			# use the basename of the rnk file to name the list elements
			n <- sapply(res, function(x) sub(".rnk", "", basename(x$rpt$rnk)))
			names(res) <- n
		}
		cat("\n")
		res
	}
	else if( is.gsea.dir(dir) ) {
		# then there was a single GSEA dir
		res <- list()
		res$tt <- import.gsea.topTable(dir)
		res$leading.edge <- import.gsea.leadingedge(dir)
		res$rnk <- import.gsea.rnk(dir)
		res$rpt <- import.gsea.rpt(dir)

		if( edb ) {
			res$edb <- import.gsea.edb(dir)
			res$gmt <- import.gsea.gmt(res$rpt$gmt)
		}
		res
	}
	else if( is.dir(dir) && any(is.gsea.dir(dir(dir, full=TRUE))) ){
		# then there must be a dir containing multiple GSEA dirs.
		dirs <- dir(dir, full=TRUE)
		dirs <- dirs[is.gsea.dir(dirs)]
		if( length(dirs) == 0 ) {
			stop("No GSEA dirs found.\n")
		}
		else {
			cat(sprintf("Found %d GSEA sub-directories in '%s'.\n", length(dirs), dir))
		}
		res <- import.gsea(dirs, edb=edb)
		res
	}
	else {
		stop("Unsupported dir argument.\n")
	}
}
