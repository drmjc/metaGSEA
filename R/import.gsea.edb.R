#' Import a GSEA edb file.
#' The edb file is a machine readable format of a GSEA run, which the GSEA
#' desktop tool can open. 
#' This function is only really indented to be run by \code{\link{import.gsea}} with 
#' \code{edb=TRUE}, and it's a good idea to import this file if you want to
#' start manipulating, subsetting, or mergeing GSEA runs and still wish to use
#' the GSEA tool (such as the GSEA Leading Edge tool) in the future. 
#' 
#' These edb files are quite large, eg 22 MB for c2_all with 1000 permutations. Size
#' is proportional to the number of permutations.
#' 
#' The edb file is quite simple XML, utilising one node per geneset, and lots
#' of attributes but not sub-terms for each geneset. We use the \code{XML} package,
#' originally from Duncan Temple Lang.
#' 
#' @note
#' Here's some code I found useful in deciphering the XML document: 
#' \code{ 
#' edb.file <- "T47D_TreatedVsUntreated/c2_all.GseaPreranked.1236217933710/edb/results.edb"
#' doc <- xmlTreeParse(edb.file)
#' r <- xmlRoot(doc); xmlAttrs(r); xmlSize(r); xmlSize(r[[1]])
#' str(xmlAttrs(r[[1]]))
#' xmlAttrs(r[[1]])[3]
#' # modify the attributes... 
#' a <- xmlAttrs(r[[1]]) a[3] <- p(a[3], "_MJC")
#' n <- xmlNode("DTG", attrs=a) r[[1]] <- n 
#' }
#' 
#' @param x A GSEA directory.
#' 
#' @return The edb file in the form of an XMLDocument. See \code{\link[XML]{xmlTreeParse}},
#'   This object appears to be immutable.
#' 
#' @author Mark Cowley, 2009-10-12
#' 
#' @export
#' @importFrom XML xmlTreeParse
#' 
#' @seealso \code{\link[XML]{xmlTreeParse}}, \code{\link{import.gsea}},
#'   \code{\link{gsea.rename.genesets.edb}}
#' 
#' @references 
#' Omega hat XML: \url{http://www.omegahat.org/RSXML/} 
#' Quick Guide to the XML package: \url{http://www.omegahat.org/RSXML/Tour.pdf}
#' @keywords file IO
#' @examples
#' \dontrun{
#' edbT <- import.gsea.edb("T47D_TreatedVsUntreated/c2_all.GseaPreranked.1236217933710")
#' edbM <- import.gsea.edb("MCF7_TreatedVsUntreated/c2_all.GseaPreranked.1236217933720")
#' edb.merged <- merge.gsea.edb("T47D_TreatedVsUntreated/c2_all.GseaPreranked.1236217933710", "MCF7_TreatedVsUntreated/c2_all.GseaPreranked.1236217933720", c("T47D_", "MCF7_"), c("", ""))
#' export.gsea.edb(edb.merged, "./merged")
#' }
import.gsea.edb <- function(x) {
	if( is.gsea.dir(x) ) {
		rpt <- import.gsea.rpt(x)
		edb.file <- rpt$edb
		if( !file.exists(edb.file) )
			stop("Can't find the edb directory. Have you moved this GSEA directory, and not updated the rpt file?")
		return( import.gsea.edb(edb.file) )
	}
	else if( is.file(x) && file.exists(x) ) {
		edb <- xmlTreeParse(x)
		edb
	}
	else {
		stop("x must be either a gsea directory, or the path to an edb file.\n")
	}
}
